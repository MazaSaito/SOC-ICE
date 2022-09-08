;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; SOC-ICE-01 Sample Code
;;
;; To prepare and run on linux:
;;  1. cp SOC_ICE_sample.tar.gz "the top directory"
;;  2. cd "the top directory"
;;  3. tar zxvf SOC_ICE_sample.tar
;;  4. cd sample
;;  5. execute "./run_SOC_ICE_sample.sh"
;;
;; to run in IDL, ".run SOC_ICE" at the IDL prompt
;;
;; SET UPs (inside the main program)
;;  * Need to set 'dir_work0' to "the top directory"
;;  * Can change parameters in 'para_file' and 'init_file' files
;;       in 'dir_init_para' 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;Import function
;;Carbon
@./func_sample/func_Kcrit
@./func_sample/func_Kn
@./func_sample/func_tau
;;Permafrost Pattern
@./func_sample/func_alt
@./func_sample/func_calc_Tamp
@./func_sample/func_n_factor
@./func_sample/func_calc_TFDD
@./func_sample/func_judge_pattern
;;Water
@./func_sample/func_gzai1
@./func_sample/func_gzai2
@./func_sample/func_init_I0W0
@./func_sample/func_Kt_dry
@./func_sample/func_Kt
@./func_sample/func_Kf
@./func_sample/func_dt
@./func_sample/func_df
@./func_sample/func_kai_Ihb
@./func_sample/func_rho
@./func_sample/func_Wo
@./func_sample/func_phy
;;LF
@./func_sample/func_make_LFparam
@./func_sample/func_LF
;;common
@./func_sample/func_r_array
@./func_sample/func_inner_product
@./func_sample/func_search_before_idx
@./func_sample/func_write
;;;;-----------------------------------------------
;;Main program
;;-------------------------------------------------
pro SOC_ICE
  print, 'start time=', systime()
  now_time = string(systime())
  year_0nen = 1950
;;;;;;;;;;;;;;;;;;;;;;
  un = '_'
  sla = '/'

;;;;;;;  
  ;; main directory
  dir_work0 = '/work01/SOC_ICE_sample/sample/'

  ;; directory where parameter files reside
  dir_init_para = dir_work0+'./data_parameter_sample/'

;;;;;;;
  ;; directory where forcing data resides
  fdir = dir_work0+'./dir_forcing_sample/data_forcing125k_SR_IPSL_UDel/'

  ;; directory where results of SOC_ICE are stored
  out_dir = dir_work0 + './data_SOC_sample/'

  ;; directories where boundary data files reside
  odir = dir_work0+'./data_misc_sample/'
  odir_f = odir
  dir_alt = odir+ 'data_Alt/'
  dir_CO2 = odir


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;setting for execution
;;read the parameter file for the years (time steps) to read in
;;the forcing (and boundary) data
  para_file = 'parameter_years_sample.csv'
  in_file = dir_init_para+para_file
  print, 'in_file parameter=', in_file
  dcsv = $
     read_ascii(in_file, $
                count = dcnt, data_start = 1, header = hd, $
                delimiter = ',', missing_value = missing)
  time_selectls = dcsv.(0)[*]
  flg_restart = 0
  help, time_selectls

  yr_ago = 0

  method_TaPr = 'SR_IPSL_UDel'
  print, 'method_TaPr =', method_TaPr

  idx_time_125 = where(time_selectls[*] eq yr_ago, cnt_time_125)
  print, 'yr_ago=', yr_ago, idx_time_125
  str_yr_ago = strtrim(string(yr_ago, format = '(i0)'), 2)
  if (cnt_time_125 ne 1) then begin
    print, 'err time should be -125k-64. check and re-input appropriate value'
    goto, skip
  endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;setting for execution
  yr_SR0 = 1950
  
  Dir3_1 = ''
  year_now = 2014
  print, 'fdir=', fdir

  fdir1 = fdir
  print, 'fdir1=', fdir1

  yr_end = year_now-yr_SR0
  print, 'yr_now yr_end=', year_now, yr_end
  print, 'Dir3_1=', Dir3_1

  no_yr_after_0nen = year_now - year_0nen

  missing = -999
  unit_mm = 1000
  sav_name = '.sav'
  nc_name = '.nc'
  csv_name =  '.csv'
  
  filename0 = 'odata'
  filename_alt0 = 'alt'
  beg_lat_st = 25.5
  beg_lat = 50.5
  diff_beg_lat_para = beg_lat - beg_lat_st

  str_scale = 10

  no_yr = 125000L
  dlt_t = 1.

;;comon parameters
  rate_Ta = -6.49e-3 ;;cesius/m
;;parameters for calc SOC and Water and Ice
  no_lp0 = 5002
  no_lp1 = no_lp0-1
  Kt0 = 1.2
  thresh_ice_Thaw = 1.e-4
  rho_i = 917      ;;kg/m3
  rho_w = 1000     ;;kg/m3
  h_b = 3000       ;;mm
  ;;sigma = 0.55     ;;Kangeki-ritu porosity m/m

;;
;;calc SOC(Carbon) & PF pattern & Water
;;

  p_namels = ['Anaktuvuk']
  p_lonlatls = [ $
               [-150.86909, 69.27636] $
               ]
  no_place = n_elements(p_namels)
  Dir1 = 'points'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;read initial and controlling parameters
;; Initial SOC, Initial Water, Tau, Beta, Eata_SOC, Eata_ICE,
;; minimumK, alpha
  init_file = 'initial_parameter_sample.csv'
  in_file = dir_init_para+init_file
  print, 'in_file initial parameter=', in_file
  dcsv = $
     read_ascii(in_file, $
                count = dcnt, data_start = 1, header = hd, $
                delimiter = ',', missing_value = missing)
  initSOC = dcsv.(0)[0, *]
  initW = dcsv.(0)[1, *]
  initTau = dcsv.(0)[2, *]
  beta = dcsv.(0)[3, *]
  eata0_SOC = dcsv.(0)[4, *]
  eata0_ICE = dcsv.(0)[5, *]
  minK = dcsv.(0)[6, *]
  alpha = dcsv.(0)[7, *]
  str_tau = 'tau'+strtrim(string(initTau, format = '(i)'), 2)
  print, 'initial SOC=', initSOC  ;; initial value of SOC
  print, 'initial water=', initW  ;; initial value of soil moisture
  print, 'tau=', initTau          ;; itial (base) value of Tau
  print, 'beta=', beta            ;; ice-water phase change efficiency
  print, 'eata0_SOC=', eata0_SOC  ;; Survival rate of SOC after glacial ablation
  print, 'eata0_ICE=', eata0_ICE  ;; Survival rate of ICE after glacial ablation
  print, 'minK=', minK            ;; minimum value of kappa (desomposition rate)
  print, 'alpha=', alpha          ;; mobility of liquid water under ice sheets

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;normal start or restart

  Dir2 = 'master'
  print, 'Dir2=', Dir2

;; READING boundary data IN
;restore lonlat
  filename = 'lonlat_GT125k.sav'
  restore, file = filepath(filename, root_dir = odir_f)
  print, 'lonlat file=', odir_f, filename
  lon_GT0 = temporary(lon_GT)
  lat_GT0 = temporary(lat_GT)
  lonF_GT0 = temporary(lonF_GT)
  latF_GT0 = temporary(latF_GT)
  grid_lon_GT0 = lon_GT0[1]-lon_GT0[0]
  grid_lat_GT0 = lat_GT0[1]-lat_GT0[0]
  adj_grid_lon_GT0 = grid_lon_GT0 /2.
  adj_grid_lat_GT0 = grid_lat_GT0 /2.
  idx_start = where(LAT_GT0[*] eq beg_lat, cnt_start)
  help, lon_GT0, lat_GT0

;restore ratio
  filename = 'ratio_GT125k_sample.sav'
  restore, file = filepath(filename, root_dir = odir_f)
  print, 'r ratio file=', odir_f, filename
  help, r_land0, r_ocean0, r_icecover0

;restore distance
  filename = 'dist_GT125k_sample.sav'
  restore, file = filepath(filename, root_dir = odir_f)
  print, 'dist file=', odir_f, filename
  help, dist0, oaff0, mapid0

;restore time list for altitude
  filename = 'total_area_sample.sav'
  restore, file = filepath(filename, root_dir = dir_alt)
  print, 'area file=', dir_alt, filename
  help, time_altls
  
;restore CO2 [ppmv]
  filename = 'CO2_edc_gtmip_RCP85.sav'
  restore, file = filepath(filename, root_dir = dir_CO2)
  print, 'CO2 file=', dir_CO2, filename
  CO2_00 = temporary(dat_CO2)
  idx_yr_CO2 = where(CO2_00[*, 1] eq year_0nen)
  idx_beg_yr_CO2 = idx_yr_CO2-no_yr
  idx_end_yr_CO2 = idx_yr_CO2+no_yr_after_0nen
  CO2_time = CO2_00[idx_beg_yr_CO2:idx_end_yr_CO2, 0]
  CO2_dat = CO2_00[idx_beg_yr_CO2:idx_end_yr_CO2, 2]
  no_time_CO2 = n_elements(CO2_time)
  no_yr_rev = idx_end_yr_CO2-idx_beg_yr_CO2
  no_yr_plus = no_yr_rev + 1
  help, CO2_time, CO2_dat, CO2_00

;restore soil
  filename = 'soil_GLDAS_convrt2GT_sample.sav'
  restore, file = filepath(filename, root_dir = odir)
  print, 'soil file=', odir, filename

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  idx_i0_beg = 0L
  idx_time_CO2_125 = 0
  idx_time_125000 = [missing]
  no_yr_plus2 = no_yr_plus

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  cnt = 0L
  idx_zero = 0L
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  foreach each_p, p_namels, idx_p do begin 
    if (p_lonlatls[0, idx_p] < 0) then p_lon = p_lonlatls[0, idx_p]+360. $
    else p_lon = p_lonlatls[0, idx_p]
    
    beg_lon_GT0 = p_lon - adj_grid_lon_GT0
    end_lon_GT0 = p_lon + adj_grid_lon_GT0
    beg_lat_GT0 = p_lonlatls[1, idx_p] - adj_grid_lat_GT0
    end_lat_GT0 = p_lonlatls[1, idx_p] + adj_grid_lat_GT0
    idx_lon_GT0 = where(beg_lon_GT0 le lon_GT0[*] and lon_GT0[*] lt end_lon_GT0, cnt_lonGT0)
    idx_lat_GT0 = where(beg_lat_GT0 le lat_GT0[*] and lat_GT0[*] lt end_lat_GT0, cnt_lat_GT0)
    lon = lon_GT0[idx_lon_GT0]
    str_lon = strtrim(string(lon*str_scale, format = '(i04)'), 2)
    lat = LAT_GT0[idx_lat_GT0]
    str_lat = strtrim(string(lat*str_scale, format = '(i03)'), 2)
    
    idx_lon = idx_lon_GT0
    idx_lat = idx_lat_GT0
    idx_lat01 = idx_lat
    print, idx_lon, idx_lat01

    soil_made = soil
    print, 'soil check=', each_p, soil, soil_made

    if (idx_zero eq 0) then begin
      para_C_LAND = missing
      para_W_LAND = missing
      para_W_ICECOVER = missing
      para_I_LAND = missing
      para_K = missing
    endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;restore Ta Pr
    filename0 = 'odata'
    filename1 = filename0+un+str_lat+un+str_lon
    filename = filename1+sav_name
    fname_TaPr = fdir1 + filename
    print, 'fname_TaPr=', fname_TaPr
    if (FILE_TEST(fname_TaPr) eq 0) then begin
      print, 'file TaPr not found =', fname_TaPr
      goto, skip_culc
    endif
    restore, fname_TaPr

    print, 'restore file hist=', fname_TaPr
    
    time_standard = temporary(TIME)
    time_ice = temporary(TIME_ICE)
    icethick = temporary(ICETHICK)
    dmy = temporary(ALT)
    dmy = temporary(ICECOVER)
    dmy = temporary(LAND)
    ta = temporary(TA)
    pr = temporary(PR)
    time_TaPr = time_standard

    time_ice = -1.*time_ice
    no_time_ice = n_elements(time_ice)

    idx_time_ice0 = func_search_before_idx(0, time_TaPr[0], time_ice)
    idx_restart0 = 0
    idx_restart1 = n_elements(time_TaPr)-1
    
    idx_restart_beg = idx_restart0
    idx_restart_end = idx_restart1
    no_time_TaPr = n_elements(time_TaPr[idx_restart_beg:idx_restart_end])

    Ta0 = temporary(ta[idx_restart_beg:idx_restart_end])
    Pr0 = temporary(pr[idx_restart_beg:idx_restart_end])

    tlen = n_elements(time_ice)
    dmy = []
    help, Ta0, Pr0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;make array
;; make array
    time125k = fltarr(no_yr_plus)+missing
    hgt_icethick = fltarr(no_yr_plus)+missing
    r_land = fltarr(no_yr_plus)+missing
    r_ocean = fltarr(no_yr_plus)+missing
    r_icecover = fltarr(no_yr_plus)+missing
    dist = fltarr(no_yr_plus)+missing
    oaff = fltarr(no_yr_plus)+missing
    mapid = fltarr(no_yr_plus)+missing
    C_land = fltarr(no_yr_plus)+missing
    C_ocean = fltarr(no_yr_plus)+missing
    C_icecover = fltarr(no_yr_plus)+missing
    soc = fltarr(no_yr_plus)+missing
    I_land = fltarr(no_yr_plus)+missing
    I_ocean = fltarr(no_yr_plus)+missing
    I_icecover = fltarr(no_yr_plus)+missing
    ICE = fltarr(no_yr_plus)+missing
    W_land = fltarr(no_yr_plus)+missing
    W_ocean = fltarr(no_yr_plus)+missing
    W_icecover = fltarr(no_yr_plus)+missing
    water = fltarr(no_yr_plus)+missing
    numPF = fltarr(no_yr_plus)+missing
    critk = fltarr(no_yr_plus)+missing
    k = fltarr(no_yr_plus)+missing
    tau = fltarr(no_yr_plus)+missing
    lf = fltarr(no_yr_plus)+missing
    CO2 = fltarr(no_yr_plus)+missing
    Ta = fltarr(no_yr_plus)+missing
    Pr = fltarr(no_yr_plus)+missing
    Tamp = fltarr(no_yr_plus)+missing
    TDD = fltarr(no_yr_plus)+missing
    FDD = fltarr(no_yr_plus)+missing
    FDD_snow = fltarr(no_yr_plus)+missing
    gzai_1 = fltarr(no_yr_plus)+missing
    gzai_2 = fltarr(no_yr_plus)+missing
    Evap = fltarr(no_yr_plus)+missing
    Qs = fltarr(no_yr_plus)+missing
    Kt_dry = fltarr(no_yr_plus)+missing
    Kt = fltarr(no_yr_plus)+missing
    Kf = fltarr(no_yr_plus)+missing
    dt = fltarr(no_yr_plus)+missing
    df = fltarr(no_yr_plus)+missing
    phy = fltarr(no_yr_plus)+missing
    rho_s = fltarr(no_yr_plus)+missing
    rho_l = fltarr(no_yr_plus)+missing
    alt = fltarr(no_yr_plus)+missing 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;restore alt
    filename_alt1 = filename_alt0+un+str_lat+un+str_lon
    filename = filename_alt1+sav_name
    print, filename
    restore, file = filepath(filename, root_dir = dir_alt)
    altls = temporary(altdiff)
    help, altls

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;spinup
;;make array for spinup
    time_spnnp = fltarr(no_lp1)+missing
    SOC_spnnp = fltarr(no_lp1)+missing
    Water_spnnp = fltarr(no_lp1)+missing
    ICE_spnnp = fltarr(no_lp1)+missing
    critk_spnnp = fltarr(no_lp1)+missing
    k_spnnp = fltarr(no_lp1)+missing
    Kt_dry_spnnp = fltarr(no_lp1)+missing
    Kt_spnnp = fltarr(no_lp1)+missing
    Kf_spnnp = fltarr(no_lp1)+missing
    dt_spnnp = fltarr(no_lp1)+missing
    df_spnnp = fltarr(no_lp1)+missing
    phy_spnnp = fltarr(no_lp1)+missing
    rho_s_spnnp = fltarr(no_lp1)+missing
    rho_l_spnnp = fltarr(no_lp1)+missing
    C_land_spnnp = fltarr(no_lp1)+missing
    C_ocean_spnnp = fltarr(no_lp1)+missing
    C_icecover_spnnp = fltarr(no_lp1)+missing
    W_land_spnnp = fltarr(no_lp1)+missing
    W_ocean_spnnp = fltarr(no_lp1)+missing
    W_icecover_spnnp = fltarr(no_lp1)+missing
    I_land_spnnp = fltarr(no_lp1)+missing
    I_ocean_spnnp = fltarr(no_lp1)+missing
    I_icecover_spnnp = fltarr(no_lp1)+missing
    gzai_1_spnnp = fltarr(no_lp1)+missing
    gzai_2_spnnp = fltarr(no_lp1)+missing
    Evap_spnnp = fltarr(no_lp1)+missing
    Qs_spnnp = fltarr(no_lp1)+missing

;;soc water initial parameter
    SOC0 = initSOC
    W00 = initW

    no_yr_tlen0 = tlen
    yr_ago_changed = -1210
    
    if (idx_p eq 0) then begin
      no_yr_all = no_yr_after_0nen + no_yr +1
      no_yr_tlen = no_yr_tlen0 - 1
    endif

;;if all data is missing, goto culc part
    idx_TaPr = where(Ta0[*] eq missing or Pr0[*] eq missing, cnt_TaPr)
    if (cnt_Tapr eq no_time_TaPr) then begin
      no_I0W0 = 2
      I0W0 = intarr(no_I0W0) + missing
      Ta = fltarr(no_yr_all)+missing
      Pr = fltarr(no_yr_all)+missing
      alt = fltarr(no_yr_all)+missing

      print, 'err Ta or PR == all data missing'
      goto, skip_culc
    endif

;;loop i0 tlen
    idx_year = 0L
    i_cnt0 = 0L
    i_cnt_all = 0L

    i0_1 = idx_i0_beg

    for i0 = 0, no_time_TaPr-1 do begin
      if i0 eq 0 then no_lp = no_lp0 else no_lp = 1

;;SOC

;;r land ocean ice n
      idx_ice = i0+idx_time_ice0
      if (idx_ice le no_time_ice-1) then begin
        r_land1 = r_land0[idx_ice]
        r_ocean1 = r_ocean0[idx_ice]
        r_icecover1 = r_icecover0[idx_ice]

        hgt_icethick1 = icethick[i0_1]
  
        dist1 = dist0[idx_ice]
        oaff1 = oaff0[idx_ice]
        mapid1 = mapid0[idx_ice]
        alt1 = func_alt(time_ice[i0_1], time_altls, altls)

      endif else begin
        r_land1 = r_land0[no_time_ice-1]
        r_ocean1 = r_ocean0[no_time_ice-1]
        r_icecover1 = r_icecover0[no_time_ice-1]
        hgt_icethick1 = icethick[no_time_ice-1]
  
        dist1 = dist0[no_time_ice-1]
        oaff1 = oaff0[no_time_ice-1]
        mapid1 = mapid0[no_time_ice-1]
        alt1 = func_alt(time_ice[no_time_ice-1], time_altls, altls)
      endelse ;;idx_ice <= no_time_ice-1
      r_land1 = r_land1[0]
      r_ocean1 = r_ocean1[0]
      r_icecover1 = r_icecover1[0]
      hgt_icethick1 = hgt_icethick1[0] 
      dist1 = dist1[0]
      oaff1 = oaff1[0]
      mapid1 = mapid1[0]
      alt1 = alt1[0]

;;r_land ocean icecover mixed ratio of area
      r_land1_mix = r_land1*(1.-r_icecover1)
      r_ocean1_mix = r_ocean1
      r_icecover1_mix = r_land1*r_icecover1

;;r land ocean ice n-1
      if (idx_ice eq 0) then begin
        r_land1_1 = r_land1
        r_ocean1_1 = r_ocean1
        r_icecover1_1 = r_icecover1
      endif else begin
        if (idx_ice le no_time_ice-1) then begin
          r_land1_1 = r_land0[idx_ice-1]
          r_ocean1_1 = r_ocean0[idx_ice-1]
          r_icecover1_1 = r_icecover0[idx_ice-1]
        endif else begin
          r_land1_1 = r_land0[no_time_ice-1]
          r_ocean1_1 = r_ocean0[no_time_ice-1]
          r_icecover1_1 = r_icecover0[no_time_ice-1]
        endelse ;;idx_ice <= no_time_ice-1
      endelse
      r_land1_1 = r_land1_1[0]
      r_ocean1_1 = r_ocean1_1[0]
      r_icecover1_1 = r_icecover1_1[0]

      r_area = func_r_array(r_land1, r_land1_1, $
                            r_ocean1, r_ocean1_1, $
                            r_icecover1, r_icecover1_1)

;;altitude
      if (alt1 gt 0) then begin
        decrease_Ta = alt1 * rate_Ta
      endif else begin
        decrease_Ta = 0.
      endelse ;;alt1 > 0

;;Ta Pr
      ta1 = Ta0[i0] + decrease_Ta
      pr1 = Pr0[i0]

      if (ta0[i0] eq missing or pr0[i0] eq missing) then flg_tapr = 1 $
      else flg_tapr = 0

;;LF
      idx_CO2_time = idx_time_CO2_125+idx_year
      CO2_1 = CO2_dat[idx_CO2_time]

      LFparam = func_make_LFparam(CO2_1)
      lf0 = func_LF(ta1, pr1, LFparam)

      for idx_lp = 0, no_lp-1 do begin ;;pre loop
      
;;;Water 
;;calc Tamplitude input:Taverage,sistance,longitude,
;;latitudeoutput:Tamp
        if (idx_ice le no_time_ice-1) then begin
          Tamp1 = func_calc_Tamp(ta1, dist1, lon, lat, time_ice[i0_1])
        endif else begin
          Tamp1 = func_calc_Tamp(ta1, dist1, lon, lat, time_ice[no_time_ice-1])
        endelse ;;idx_ice <= no_time_ice-1
     
;;;calc TDD and FDD from Tav and Tamp
        n_factor = func_n_factor(pr1, missing)
        TFDD = func_calc_TFDD(ta1, Tamp1) ;;input:Taverage,Tamplitude output:[TDD,FDD]
        TDD1 = TFDD[0]
        FDD1 = TFDD[1]
        FDD1_snow = n_factor * FDD1

;;;define permafrost pattern
        num_pf = func_judge_pattern(TDD1, FDD1)
        case num_pf of
          1:begin
            str_pf = 'CP'
          end
          2:begin
            str_pf = 'EP'
          end
          3:begin
            str_pf = 'SF'
          end
          4:begin
            str_pf = 'IF'
          end
          5:begin
            str_pf = 'NF'
          end
        endcase

;; define initial I0 W0
        if i0 eq 0 then begin
          if (idx_lp eq 0) then begin
            I0W0 = func_init_I0W0(W00, num_pf)
            no_I0W0 = n_elements(I0W0)
            I00 = I0W0[0]
            W0 = I0W0[1]
            Wice0 = W0
          endif ;;idx_lp == 0 
        endif   ;;i0 == 0

;;;calc gzai
        gzai_00 = func_gzai1(num_pf, soil_made, ta1, missing)
        Evap1 = gzai_00[0]
        Qs1 = gzai_00[1]
        Qsb1 = gzai_00[2]
        gzai1 = gzai_00[3]

        if (idx_ice le no_time_ice-1) then begin
          if(time_ice[i0_1] le -20100) then deltaT = 100. $
          else if(time_ice[i0_1] le -10050) then deltaT = 50. $
          else if(time_ice[i0_1] le yr_ago_changed) then deltaT = 10. $
          else deltaT = 1.
        endif else begin
          deltaT = 1.
        endelse ;;idx_ice <= no_time_ice-1


        if (idx_lp eq no_lp-1) then begin
          if (i0 eq tlen-1) then begin
            idx1 = 0
            idx2 = 0
          endif else begin
            idx1 = 0
            idx2 = fix(deltaT-1)
            
            tmp_SOC = fltarr(deltaT)+missing
            tmp_ICE = fltarr(deltaT)+missing
            tmp_W = fltarr(deltaT)+missing
            
            tmp_C_land = fltarr(deltaT)+missing
            tmp_C_ocean = fltarr(deltaT)+missing
            tmp_C_icecover = fltarr(deltaT)+missing
            
            tmp_W_land = fltarr(deltaT)+missing
            tmp_W_ocean = fltarr(deltaT)+missing
            tmp_W_icecover = fltarr(deltaT)+missing
          
            tmp_I_land = fltarr(deltaT)+missing
            tmp_I_ocean = fltarr(deltaT)+missing
            tmp_I_icecover = fltarr(deltaT)+missing
          endelse ;;i0 == tlen-1
        endif else begin
          idx1 = 0
          idx2 = 0
        endelse ;;idx_lp == no_lp-1


        for i1 = idx1, idx2 do begin

  ;;Water under ocean
          if(idx_year eq 0) then begin
            Wo_00 = func_Wo(h_b, soil_made)
        
            sigma_hb = Wo_00[0]
            Wo1 = sigma_hb
            max_W = Wo_00[1]
            min_W = Wo_00[2]
          endif ;;idx_year == 0

;;;SOC calculation
      ;;tau
          W_rate =  W0 / max_W

          tau1 = func_tau(num_pf, initTau, W_rate, missing)

      ;;critical K K
          critK1 =  func_Kcrit(SOC0, lf0, missing)
          if(flg_restart eq 0 and $
             i0 eq 0 and i1 eq 0 and idx_lp eq 0) then kk00 = critK1 $
          else kk00 = kk0
          kk0 = func_Kn(critK1, kk00, tau1)
          if (kk0 lt minK) then kk0 = mink

      ;;SOC icecover filter
          dcmpc0 = kk0 * SOC0
          dcmpc1 = lf0 - dcmpc0
          SOC1 = SOC0 + dlt_t * dcmpc1
     
;;;Water calculation
     ;;limit Ice contribution for Kf-evaluation
          I_hb = I00/h_b
          kai_Ihb1 = func_kai_Ihb(I00, h_b, I_hb, sigma_hb)

      ;;phy
          Kt_dry1 = func_Kt_dry(SOC0, h_b, sigma_hb)
          Kt1 = func_Kt(SOC0, h_b, sigma_hb, W0, Kt_dry1)
          Kf1 = func_Kf(Kt1, kai_Ihb1)
          rho = func_rho(rho_i, kai_Ihb1, rho_w, W0, h_b)
          rho_s1 = rho[0]
          rho_l1 = rho[1]
          rho1 = rho_s1+rho_l1
          dt1 = func_dt(Kf1, TDD1, rho1)
          df1 = func_df(Kt1, FDD1_snow, rho1)
          phy1 = func_phy(df1, dt1, beta, W0, I00, h_b, num_pf, missing)
          phy1_ice = func_phy(df1, dt1, beta, Wice0, I00, h_b, num_pf, missing)

     ;;Water Wn current from landice melt
          W_i_lio = 0.
          if (idx_year gt 0) then begin
            if (r_icecover[idx_year-1] gt thresh_ice_thaw) then begin
              if (r_icecover[idx_year] le thresh_ice_thaw) then begin
                idx_Ice_thaw = i0+idx_time_ice0
                if (idx_Ice_thaw ge no_time_ice-1) then idx_Ice_thaw = no_time_ice-1
                icethick_melt = icethick[idx_Ice_thaw]
                W_i_lio = icethick_melt * unit_mm / deltaT
              endif ;;r_icecover[idx_year]
            endif   ;;r_icecover[idx_year-1]
          endif     ;;r_ice_year

          if (i_cnt_all eq 0) then Wice0 = W0
          W_i_lio1 = W_i_lio * r_icecover[idx_year-1]
          gzai2 = func_gzai2(soil_made, missing)

      ;;Water land
          W1 = ((1. - gzai1) * pr1) + ((1. - gzai2) * W0) - phy1
          if (W1 gt max_W) then begin
            W1 = max_W
          endif else if (W1 lt 0.) then begin
            print, 'err w1 <0) =', i0, i1, W1
            W1 = 0.
          endif else W1 = W1

      ;;Water under icecover
          Wice1 = (alpha * (1. - gzai1) * W_i_lio1) + ((1. - gzai2) * Wice0) - phy1_ice
          if (Wice1 gt max_W) then begin
            Wice1 = max_W
          endif else if (Wice1 lt 0.) then begin
            print, 'err wice1 <0) =', i0, i1, Wice1
            Wice1 = 0.
          endif else Wice1 = Wice1

      ;;W_i_lio
          if (W_i_lio eq 0) then begin
            eata_SOC = 1.0
            eata_ICE = 1.0
          endif else if (W_i_lio eq 1) then begin
            eata_SOC = eata0_SOC
            eata_ICE = eata0_ICE
          endif else begin
            eata_SOC = 1.0
            eata_ICE = 1.0
          endelse
          if (Wice1 lt 0.) then begin
            print, 'Wice1=', Wice1, W_i_lio1, W_i_lio, Wice0, phy1_ice
          endif

;;C land ocean ice 
      ;;land
          C_loi = [SOC1, SOC0, eata_SOC*SOC0]
          tmp_C_land1 = func_inner_product(r_area[*, 0], C_loi)
      ;;ocean
          C_loi = [SOC0, SOC0, eata_SOC*SOC0]
          tmp_C_ocean1 = func_inner_product(r_area[*, 1], C_loi)
      ;;ice
          C_loi = [SOC0, SOC0, SOC1]
          tmp_C_icecover1 = func_inner_product(r_area[*, 2], C_loi)

          if (tmp_C_land1 lt 0.0) then begin
            print, 'err tmp_C_land =', tmp_C_land1, r_area, C_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif
          if (tmp_C_ocean1 lt 0.0) then begin
            print, 'err tmp_C_ocean =', tmp_C_ocean1, r_area, C_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif
          if (tmp_C_icecover1 lt 0.0) then begin
            print, 'err tmp_C_icecover =', tmp_C_icecover1, r_area, C_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif

          if (tmp_C_land1 lt 0.0) then tmp_C_land1 = 0.
          if (tmp_C_ocean1 lt 0.0) then tmp_C_ocean1 = 0.
          if (tmp_C_icecover1 lt 0.0) then tmp_C_icecover1 = 0.

          tmp_SOC1 = $
             (r_land1_mix*tmp_C_land1)+ $
             (r_ocean1_mix*tmp_C_ocean1)+ $
             (r_icecover1_mix*tmp_C_icecover1)
        
          if (idx_lp eq no_lp-1) then begin
    ;;make array
            tmp_SOC[i1] = tmp_SOC1
            tmp_C_land[i1] = tmp_C_land1
            tmp_C_ocean[i1] = tmp_C_ocean1
            tmp_C_icecover[i1] = tmp_C_icecover1
            
            if (SOC1 lt 0.) then SOC1 =  0.
            if (flg_tapr eq 0) then begin
              Ta[idx_year] = ta1
              Pr[idx_year] = pr1
              critK[idx_year] = critK1
              k[idx_year] = kk0
              lf[idx_year] = lf0
              tau[idx_year] = tau1
            endif ;;flg_tapr=0

            r_land[idx_year] = r_land1
            r_ocean[idx_year] = r_ocean1
            r_icecover[idx_year] = r_icecover1

            hgt_icethick[idx_year] = hgt_icethick1
            dist[idx_year] = dist1        
            oaff[idx_year] = oaff1
            mapid[idx_year] = mapid1
            alt[idx_year] = alt1

          endif ;;idx_lp


      ;;ICE
          I11 = I00 + phy1
          if (I11 lt 0.) then I11 = 0.

;;I land ocean ice 
      ;;land
          I_loi = [I11, I00, eata_ICE*I00]
          tmp_I_land1 = func_inner_product(r_area[*, 0], I_loi)
      ;;ocean
          I_loi = [I00, I00, eata_ICE*I00]
          tmp_I_ocean1 = func_inner_product(r_area[*, 1], I_loi)
      ;;ice
          I_loi = [I00, I00, I11]
          tmp_I_icecover1 = func_inner_product(r_area[*, 2], I_loi)

          if (tmp_I_land1 lt 0.0) then begin
            print, 'err tmp_I_land =', tmp_I_land1, r_area, I_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif
          if (tmp_I_ocean1 lt 0.0) then begin
            print, 'err tmp_I_ocean =', tmp_I_ocean1, r_area, I_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif
          if (tmp_I_icecover1 lt 0.0) then begin
            print, 'err tmp_I_icecover =', tmp_I_icecover1, r_area, I_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif

          if (tmp_I_land1 lt 0.0) then tmp_I_land1 = 0.
          if (tmp_I_ocean1 lt 0.0) then tmp_I_ocean1 = 0.
          if (tmp_I_icecover1 lt 0.0) then tmp_I_icecover1 = 0.


      ;;Ice total
          tmp_I11 = $
             r_land1_mix * tmp_I_land1 + $
             r_ocean1_mix * tmp_I_ocean1 + $
             r_icecover1_mix * tmp_I_icecover1

          if (idx_lp eq no_lp-1) then begin
            tmp_I_land[i1] = tmp_I_land1
            tmp_I_ocean[i1] = tmp_I_ocean1
            tmp_I_icecover[i1] = tmp_I_icecover1
            tmp_ICE[i1] = tmp_I11
          endif ;;idx_lp

      ;;Water total
;;I land ocean ice 
      ;;land
          W_loi = [W1, Wo1, W0]
          tmp_W_land1 = func_inner_product(r_area[*, 0], W_loi)
      ;;ocean
          W_loi = [W0, Wo1, Wice1]
          tmp_W_ocean1 = func_inner_product(r_area[*, 1], W_loi)
      ;;ice
          W_loi = [W0, Wo1, Wice1]
          tmp_W_icecover1 = func_inner_product(r_area[*, 2], W_loi)

          if (tmp_W_land1 lt 0.0) then begin
            print, 'err tmp_W_land =', tmp_W_land1, Wo1, r_area, W_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif
          if (tmp_W_ocean1 lt 0.0) then begin
            print, 'err tmp_W_ocean =', tmp_W_ocean1, Wo1, r_area, W_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif
          if (tmp_W_icecover1 lt 0.0) then begin
            print, 'err tmp_W_icecover =', tmp_W_icecover1, Wo1, r_area, W_loi
            print, 'r loi =', r_land1, r_land1_1, $
                   r_ocean1, r_ocean1_1, r_icecover1, r_icecover1_1
          endif

          if (tmp_W_land1 lt 0.0) then tmp_W_land1 = 0.
          if (tmp_W_ocean1 lt 0.0) then tmp_W_ocean1 = 0.
          if (tmp_W_icecover1 lt 0.0) then tmp_W_icecover1 = 0.


          tmp_W1 = $
             r_land1_mix * tmp_W_land1 + $
             r_ocean1_mix * tmp_W_ocean1 + $
             r_icecover1_mix * tmp_W_icecover1

          if (idx_lp eq no_lp-1) then begin
            tmp_W_land[i1] = tmp_W_land1
            tmp_W_ocean[i1] = tmp_W_ocean1
            tmp_W_icecover[i1] = tmp_W_icecover1
            tmp_W[i1] = tmp_W1
          endif ;;idx_lp

          SOC0 = tmp_C_land1
          I00 = tmp_I_land1
          W0 = tmp_W_land1
          Wice0 = tmp_W_icecover1

          if (i0 eq 0 and idx_lp ne no_lp-1) then begin
            idx_lp1 = idx_lp

            time_spnnp[idx_lp1] = no_lp1 - idx_lp1
            SOC_spnnp[idx_lp1] = tmp_SOC1
            Water_spnnp[idx_lp1] = tmp_W1
            ICE_spnnp[idx_lp1] = tmp_I11
            critK_spnnp[idx_lp1] = critK1
            k_spnnp[idx_lp1] = kk0
            Kt_dry_spnnp[idx_lp1] = Kt_dry1
            Kt_spnnp[idx_lp1] = Kt1
            Kf_spnnp[idx_lp1] = Kf1
            dt_spnnp[idx_lp1] = dt1
            df_spnnp[idx_lp1] = df1
            phy_spnnp[idx_lp1] = phy1
            rho_s_spnnp[idx_lp1] = rho_s1
            rho_l_spnnp[idx_lp1] = rho_l1
            C_land_spnnp[idx_lp1] = tmp_C_land1
            C_ocean_spnnp[idx_lp1] = tmp_C_ocean1
            C_icecover_spnnp[idx_lp1] = tmp_C_icecover1
            W_land_spnnp[idx_lp1] = tmp_W_land1
            W_ocean_spnnp[idx_lp1] = tmp_W_ocean1
            W_icecover_spnnp[idx_lp1] = tmp_W_icecover1
            I_land_spnnp[idx_lp1] = tmp_I_land1
            I_ocean_spnnp[idx_lp1] = tmp_I_ocean1
            I_icecover_spnnp[idx_lp1] = tmp_I_icecover1
            gzai_1_spnnp[idx_lp1] = gzai1
            gzai_2_spnnp[idx_lp1] = gzai2
            Evap_spnnp[idx_lp1] = Evap1
            Qs_spnnp[idx_lp1] = Qs1
          endif  ;; idx_lp

          if (idx_lp eq no_lp-1) then begin
    ;;make array
            time125k[idx_year] = 1L*time_ice[0]+idx_year
            if (flg_tapr eq 0) then begin
              numPF[idx_year] = num_pf
              Tamp[idx_year] = Tamp1
              TDD[idx_year] = TDD1
              FDD[idx_year] = FDD1
              FDD_snow[idx_year] = FDD1_snow
              gzai_1[idx_year] = gzai1
              gzai_2[idx_year] = gzai2
              Evap[idx_year] = Evap1
              Qs[idx_year] = Qs1
              Kt_dry[idx_year] = Kt_dry1
              Kt[idx_year] = Kt1
              Kf[idx_year] = Kf1
              dt[idx_year] = dt1
              df[idx_year] = df1
              phy[idx_year] = phy1
              rho_s[idx_year] = rho_s1
              rho_l[idx_year] = rho_l1
            endif ;;flg_tapr=0

            CO2[idx_year] = CO2_1
            idx_year = idx_year+1
          endif ;;idx_lp
      
          i_cnt_all = i_cnt_all+1
        endfor ;;i1
      endfor   ;;idx_lp for initial values

      if(i_cnt0 eq 0) then begin
        idx_beg = 0L
        idx_end = long((idx_beg+deltaT)-1)
      endif else begin
        idx_beg = idx_end+1
        idx_end = long((idx_beg+deltaT)-1)
      endelse

      if (flg_tapr eq 0) then begin
        C_land[idx_beg:idx_end] = temporary(tmp_C_land[*])
        C_ocean[idx_beg:idx_end] = temporary(tmp_C_ocean[*])
        C_icecover[idx_beg:idx_end] = temporary(tmp_C_icecover[*])
        soc[idx_beg:idx_end] = temporary(tmp_SOC[*])
       
        I_land[idx_beg:idx_end] = temporary(tmp_I_land[*])
        I_ocean[idx_beg:idx_end] = temporary(tmp_I_ocean[*])
        I_icecover[idx_beg:idx_end] = temporary(tmp_I_icecover[*])
        ICE[idx_beg:idx_end] = temporary(tmp_ICE[*])
       
        W_land[idx_beg:idx_end] = temporary(tmp_W_land[*])
        W_ocean[idx_beg:idx_end] = temporary(tmp_W_ocean[*])
        W_icecover[idx_beg:idx_end] = temporary(tmp_W_icecover[*])
        water[idx_beg:idx_end] = temporary(tmp_W[*])
      endif   ;;flg_tapr=0

      i_cnt0 = i_cnt0+1
      if (i0_1 ge no_time_ice-1) then i0_1 = no_time_ice-1 $
      else i0_1 = i0_1+1
    endfor                      ;tlen, i0

;;get value as selected year
    if (idx_zero eq 0) then begin
      idx_time_125000 = where(TIME125K[*] eq yr_ago, cnt_time_125000)
      print, 'idx_time 125000= ', idx_time_125000, yr_ago
      idx_time_1yr_before = idx_time_125000-1
      no_yr_plus2 = (no_yr_rev - idx_time_125000) + 1

      para_C_LAND = missing
      para_W_LAND = missing
      para_W_ICECOVER = missing
      para_I_LAND = missing
      para_K = missing
    endif

    if (cnt_time_125000 ne 1) then begin
      print, 'err. index make para data', idx_time_125000
      goto, skip
    endif
    para_C_LAND = C_LAND[idx_time_125000]
    para_W_LAND = W_LAND[idx_time_125000]
    para_W_ICECOVER = W_ICECOVER[idx_time_125000]
    para_I_LAND = I_LAND[idx_time_125000]
    para_K = K[idx_time_125000]

;;;write sav file
    fnameout = 'odata_f'+un+str_lat+un+str_lon
    dir_out0 = out_dir+Dir1+sla+method_TaPr+sla+Dir3_1+str_tau
    file_mkdir, dir_out0
    dir_out = dir_out0 + sla
    fnameout_sav = dir_out+fnameout+sav_name
    print, 'output sav =', fnameout_sav

    file_out_nc = dir_out+fnameout+nc_name
    file_out_csv = dir_out+fnameout+csv_name
    file_out_spnnp_csv = dir_out+fnameout+'_spnnp'+csv_name

;;goto, skip

;;normal start

    readme = 'time 125k-64, normal start since -125k spinnup 5000'

    func_write, $
       missing, flg_restart, readme, $
       tlen, no_yr_all, no_I0W0, no_lp1, $
       fnameout_sav, file_out_nc, file_out_csv, file_out_spnnp_csv, $
       time_ice, time125k, $
       hgt_icethick, r_land, r_ocean, r_icecover, $
       dist, oaff, mapid, $
       Ta, Pr, critK, k, tau, CO2, lf, soc, $
       Tamp, TDD, FDD, FDD_snow, numPF, $
       gzai_1, gzai_2, Evap, Qs, $
       Kt_dry, Kt, Kf, rho_s, rho_l, dt, df, phy, $
       ICE, water, I0W0, $
       C_land, C_Ocean, C_icecover, $
       W_land, W_Ocean, W_icecover, $
       I_land, I_Ocean, I_icecover, $
       alt, $
       SOC_spnnp, Water_spnnp, ICE_spnnp, $
       critK_spnnp, k_spnnp, $
       Kt_dry_spnnp, Kt_spnnp, Kf_spnnp, $
       dt_spnnp, df_spnnp, phy_spnnp, $
       rho_s_spnnp, rho_l_spnnp, $
       C_land_spnnp, C_ocean_spnnp, C_icecover_spnnp, $
       W_land_spnnp, W_ocean_spnnp, W_icecover_spnnp, $
       I_land_spnnp, I_ocean_spnnp, I_icecover_spnnp, $
       gzai_1_spnnp, gzai_2_spnnp, Evap_spnnp, Qs_spnnp, $
       time_spnnp

    idx_zero = idx_zero+1

skip_culc:
    
    if(cnt mod 200 eq 0) then print, 'cnt=', cnt
    cnt = cnt+1
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;goto, skip
  endforeach  ;;each_p
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


skip:
  print, 'end=', systime()
  print, 'Process Finished'
  return
end 
