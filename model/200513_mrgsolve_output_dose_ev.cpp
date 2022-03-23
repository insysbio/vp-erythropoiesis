$PROB
# Model: `mrg1`
  - Title: 
  - Notes: 
  - Source: Generated automatically from platform with qs3p-js 0.3.0

# Demo
```{r,echo=TRUE}
  ev(amt=10) %>% mrgsim %>% plot
```

$SET end=120, delta=0.1, hmax=0.01, hmin=0, rtol=1e-3, atol=1e-6

$PARAM @annotated
// @Const ''
time_inj : 0.01 : (h)
// @Const ''
IU_per_pM : 5.25 : (UL/pmole)
// @Const ''
Mr_epo : 35000 : (g/mole)
// @Const ''
coef_epo_sc : 600 : (UL/UL)
// @Const ''
h_epo_sc : 1.5 : (UL)
// @Const ''
k_epo_sc_lym : 0.06016045 : (1/h)
// @Const ''
k_epo_lym_pl : 0.05427484 : (1/h)
// @Const ''
k_epo_sc_pl : 0.005577687 : (1/h)
// @Const ''
k_epo_bm_pl : 10 : (1/h)
// @Const ''
k_epo_per_pl : 0.5733347 : (1/h)
// @Const ''
k_epo_pl_bm : 2.256684 : (1/h)
// @Const ''
k_epo_pl_per : 0.06634371 : (1/h)
// @Const ''
k_dis_epoR : 30 : (1/h)
// @Const ''
Kd_epoR : 45.18755 : (pM)
// @Const ''
k_int_epoR : 0.04704679 : (1/h)
// @Const ''
k_deg_epoR : 0.01018734 : (1/h)
// @Const ''
k_syn_epoR_cfu_e : 12.4095 : (item/h/cell)
// @Const ''
k_ass_epoR_pl : 0.01018734 : (1/pM/h)
// @Const ''
k_dis_epoR_pl : 0.01018734 : (1/h)
// @Const ''
k_syn_epoR_pl : 0.01018734 : (pM/h)
// @Const ''
k_deg_epoR_pl : 0.01018734 : (1/h)
// @Const ''
k_syn_hif1a : 1112.014 : (UL/h)
// @Const ''
k_syn_epo_mrna : 1.19972 : (1/h)
// @Const ''
IC50_roxa : 1.7 : (uM)
// @Const ''
HIF1a_base : 351 : (UL)
// @Const ''
gamma_hif1a : 1.794074 : (UL)
// @Const ''
gamma_hb : 3.246693 : (UL)
// @Const ''
kbase_syn_epo_pl : 0.4079528 : (pM/h)
// @Const ''
k_deg_epo_pl : 0.2063347 : (1/h)
// @Const ''
Vmax_deg_epo_pl : 34.48391 : (pM/h)
// @Const ''
Km_deg_epo_pl : 431.8562 : (pM)
// @Const 'Molecular weight of Roxadustat'
Mr_Roxa : 352.346 : (g/mole)
// @Const ''
F_Roxa : 0.4466476 : (UL)
// @Const ''
Dose_Roxa_PO_mg : 0 : (mg)
// @Const ''
k_deg_roxa_pl : 0.2726787 : (1/h)
// @Const ''
k_absorb_roxa : 0.4202667 : (1/h)
// @Const ''
k_roxa_pl_per : 0.1501216 : (1/h)
// @Const 'rate constant for transport of Roxadustat from peripheral compartment into blood plasma'
k_roxa_per_pl : 0.06076102 : (1/h)
// @Const ''
BMI : 21.7 : (UL)
// @Const ''
BW : 70 : (kg)
// @Const ''
BH : 1.8 : (m)
// @Const ''
sex : 1 : (UL)
// @Const ''
EC50_proe : 2090835 : (kcell/L)
// @Const ''
HSC_bm_base : 1428000 : (kcell/L)
// @Const ''
MPP_bm_base : 1714000 : (kcell/L)
// @Const ''
CMP_bm_base : 14283000 : (kcell/L)
// @Const ''
N_A_pmole : 602200000000 : (item/pmole)
// @Const ''
kbase_pro : 1 : (1/h)
// @Const ''
kmax_pro_cfu_e_bm : 4.935083 : (1/h)
// @Const ''
k_pro_clp_bm : 1 : (1/h)
// @Const ''
kbase_pro_bfu : 0.003567283 : (1/h)
// @Const ''
EC50_pro_epo_cfu_e : 0.023 : (item/item)
// @Const ''
EC50_dif_epo_cfu_e : 0.01 : (item/item)
// @Const ''
Emax_dif_epo : 9.089291 : (UL)
// @Const ''
coef_dif_epo : 873.7265 : (item/cell)
// @Const ''
h_dif_epo : 6.162757 : (UL)
// @Const ''
kbase_dif_bfu_le_cfu_e_bm : 0.01647444 : (1/h)
// @Const ''
kmax_dif_cfu_e_preproe_bm : 4.255028 : (1/h)
// @Const ''
kbase_dif_bm : 0.08240409 : (1/h)
// @Const ''
k_dif_hsc_mpp_bm : 0.8104074 : (1/h)
// @Const ''
k_dif_mpp_clp_bm : 0.04946078 : (1/h)
// @Const ''
k_dif_mpp_cmp_bm : 1 : (1/h)
// @Const ''
k_dif_cmp_gmp_bm : 0.3100506 : (1/h)
// @Const ''
k_dif_cmp_mep_bm : 0.001581855 : (1/h)
// @Const ''
k_dif_mpp_clp : 0.04946078 : (1/h)
// @Const ''
k_dif_mep_bfu_ee_bm : 0.01002318 : (1/h)
// @Const ''
k_dif_mep_mkc_bm : 0.01 : (1/h)
// @Const ''
k_dif_bfu_ee_bfu_le_bm : 0.02 : (1/h)
// @Const ''
k_dif_gmp_myelobl_bm : 0.02784777 : (1/h)
// @Const ''
kbase_dif_orthoe_ret_bm : 0.04272208 : (1/h)
// @Const ''
k_dif_ret_ery_pl : 0.035 : (1/h)
// @Const ''
k_migr_bm_pl : 1e-12 : (1/h)
// @Const ''
k_migr_ret_bm_pl : 0.03120928 : (1/h)
// @Const ''
kbase_apo_hsc : 0.06879673 : (1/h)
// @Const ''
kbase_apo_bm : 0.009753309 : (1/h)
// @Const ''
k_apo_clp_bm : 0.3895091 : (1/h)
// @Const ''
k_apo_gmp_bm : 0.4835958 : (1/h)
// @Const ''
k_apo_ret_bm : 0.0001 : (1/h)
// @Const ''
k_apo_mkc_bm : 1 : (1/h)
// @Const ''
k_apo_pl : 0.06 : (1/h)
// @Const ''
k_apo_ery_pl : 0.000347222222222222 : (1/h)
// @Const ''
k_apo_ret_pl : 0.0001 : (1/h)
// @Const ''
k_apo_ne_myelobl_bm : 0.007633716 : (1/h)
// @Const ''
kbase_apo_cfu_e_bm : 0.07716218 : (1/h)
// @Const ''
IC50_apo_epo_cfu_e : 0.5726052 : (UL)
// @Const ''
Imax_apo_epo : 0.0008805405 : (UL)
// @Const ''
ka_epo_bfu : 0.1 : (UL)
// @Const ''
k_scf_epo_apo_bfu : 1 : (UL)
// @Const ''
k_int_scfR : 5.4 : (1/h)
// @Const ''
k_ass_scfR : 0.017 : (1/h/pM)
// @Const ''
k_dis_scfR : 13.69215 : (1/h)
// @Const ''
k_deg_scf_pl : 0.01 : (1/h)
// @Const ''
k_deg_scf_bm : 0 : (1/h)
// @Const ''
k_deg_scfR_pl : 0.1 : (1/h)
// @Const ''
k_deg_scfR_bm : 0.7153708 : (1/h)
// @Const ''
k_prod_scf_pl : 0.5 : (pM/h)
// @Const ''
k_prod_scf_bm : 1837.919 : (pM/h)
// @Const ''
k_prod_scfR_pl : 0.5 : (pM/h)
// @Const ''
k_syn_scfR_bm : 17894050 : (item/h/kcell)
// @Const ''
k_scf_bm_pl : 0.00019904458 : (1/h)
// @Const ''
k_scf_pl_bm : 0.00003377237 : (1/h)
// @Const ''
Emax_pro_scf : 2 : (UL)
// @Const ''
Imax_apo_scf : 0.01 : (UL)
// @Const ''
Ka_scf : 0.0003 : (1/pM)
// @Const ''
PS_scfR_pl_bm : 0.01 : (L/h)
// @Const ''
k_deg_il3_pl : 5.8 : (1/h)
// @Const ''
Vm_deg_il3_pl : 432.3659 : (pM/h)
// @Const ''
Km_deg_il3_pl : 31.36738 : (pM)
// @Const ''
kcell_deg_il3_bm : 0 : (UL)
// @Const ''
kbase_deg_il3_bm : 1 : (1/h)
// @Const ''
k_prod_il3_pl : 19.32172 : (pM/h)
// @Const ''
k_prod_il3_bm : 0 : (pM/h)
// @Const ''
k_il3_bm_pl : 2.483092 : (1/h)
// @Const ''
k_il3_pl_bm : 0.588715 : (1/h)
// @Const ''
k_il3_pl_per : 15.01982 : (1/h)
// @Const ''
k_il3_per_pl : 3.57891 : (1/h)
// @Const ''
Emax_pro_il3 : 1 : (UL)
// @Const ''
EC50_pro_il3 : 20 : (pM)
// @Const ''
Imax_apo_il3 : 0.00001 : (UL)
// @Const ''
IC50_apo_il3 : 20 : (pM)
// @Const ''
epoR_bfu_le_total_exp : 80 : (item/cell)
// @Const ''
epoR_cfu_e_total_exp : 1100 : (item/cell)
// @Const ''
epoR_preproe_total_exp : 1100 : (item/cell)
// @Const ''
epoR_proe_total_exp : 1100 : (item/cell)
// @Const ''
epoR_baso_ee_total_exp : 600 : (item/cell)
// @Const ''
epoR_baso_le_total_exp : 600 : (item/cell)
// @Const ''
epoR_polye_total_exp : 320 : (item/cell)
// @Const ''
epoR_orthoe_total_exp : 80 : (item/cell)

$CMT @annotated
// @Species ''
switch_inj : amount
// @Species ''
EPO_sc : amount
// @Species ''
EPO_lym : amount
// @Species ''
EPO_mRNA_amt_ : as amount
// @Species ''
HIF1a : amount
// @Species 'Roxadustat in blood plasma'
Roxa_pl : amount
// @Species ''
Roxa_gut : amount
// @Species ''
Roxa_per : amount
// @Species ''
HSC_bm_amt_ : as amount
// @Species ''
MPP_bm_amt_ : as amount
// @Species ''
CMP_bm_amt_ : as amount
// @Species ''
CLP_bm_amt_ : as amount
// @Species ''
MEP_bm_amt_ : as amount
// @Species ''
BFU_eE_bm_amt_ : as amount
// @Species ''
BFU_lE_bm_amt_ : as amount
// @Species ''
CFU_E_bm_amt_ : as amount
// @Species ''
PreproE_bm_amt_ : as amount
// @Species ''
ProE_bm_amt_ : as amount
// @Species ''
Baso_eE_bm_amt_ : as amount
// @Species ''
Baso_lE_bm_amt_ : as amount
// @Species ''
PolyE_bm_amt_ : as amount
// @Species ''
OrthoE_bm_amt_ : as amount
// @Species ''
RET_bm_amt_ : as amount
// @Species ''
GMP_bm_amt_ : as amount
// @Species ''
MKC_bm_amt_ : as amount
// @Species ''
NE_myelobl_bm_amt_ : as amount
// @Species ''
HSC_pl_amt_ : as amount
// @Species ''
MPP_pl_amt_ : as amount
// @Species ''
CLP_pl_amt_ : as amount
// @Species ''
MEP_pl_amt_ : as amount
// @Species ''
BFU_eE_pl_amt_ : as amount
// @Species ''
BFU_lE_pl_amt_ : as amount
// @Species ''
CMP_pl_amt_ : as amount
// @Species ''
GMP_pl_amt_ : as amount
// @Species ''
RET_pl_amt_ : as amount
// @Species ''
ERY_pl_amt_ : as amount
// @Species ''
IL3_bm_amt_ : as amount
// @Species ''
IL3_pl_amt_ : as amount
// @Species ''
IL3_per_amt_ : as amount
// @Species ''
SCF_bm_amt_ : as amount
// @Species ''
SCF_pl_amt_ : as amount
// @Species ''
scfR_pl_amt_ : as amount
// @Species ''
R_SCF_bm : amount
// @Species ''
scfR_bm : amount
// @Species ''
EPO_pl_amt_ : as amount
// @Species ''
EPO_bm_amt_ : as amount
// @Species ''
EPO_per_amt_ : as amount
// @Species ''
epoR_bfu_le : amount
// @Species ''
epoR_cfu_e : amount
// @Species ''
epoR_preproe : amount
// @Species ''
epoR_proe : amount
// @Species ''
epoR_baso_ee : amount
// @Species ''
epoR_baso_le : amount
// @Species ''
epoR_polye : amount
// @Species ''
epoR_orthoe : amount
// @Species ''
R_EPO_bfu_le : amount
// @Species ''
R_EPO_cfu_e : amount
// @Species ''
R_EPO_preproe : amount
// @Species ''
R_EPO_proe : amount
// @Species ''
R_EPO_baso_ee : amount
// @Species ''
R_EPO_baso_le : amount
// @Species ''
R_EPO_polye : amount
// @Species ''
R_EPO_orthoe : amount

$GLOBAL
#define EPO_mRNA (EPO_mRNA_amt_ / Default)
#define HSC_bm (HSC_bm_amt_ / Vol_bm)
#define MPP_bm (MPP_bm_amt_ / Vol_bm)
#define CMP_bm (CMP_bm_amt_ / Vol_bm)
#define CLP_bm (CLP_bm_amt_ / Vol_bm)
#define MEP_bm (MEP_bm_amt_ / Vol_bm)
#define BFU_eE_bm (BFU_eE_bm_amt_ / Vol_bm)
#define BFU_lE_bm (BFU_lE_bm_amt_ / Vol_bm)
#define CFU_E_bm (CFU_E_bm_amt_ / Vol_bm)
#define PreproE_bm (PreproE_bm_amt_ / Vol_bm)
#define ProE_bm (ProE_bm_amt_ / Vol_bm)
#define Baso_eE_bm (Baso_eE_bm_amt_ / Vol_bm)
#define Baso_lE_bm (Baso_lE_bm_amt_ / Vol_bm)
#define PolyE_bm (PolyE_bm_amt_ / Vol_bm)
#define OrthoE_bm (OrthoE_bm_amt_ / Vol_bm)
#define RET_bm (RET_bm_amt_ / Vol_bm)
#define GMP_bm (GMP_bm_amt_ / Vol_bm)
#define MKC_bm (MKC_bm_amt_ / Vol_bm)
#define NE_myelobl_bm (NE_myelobl_bm_amt_ / Vol_bm)
#define HSC_pl (HSC_pl_amt_ / Vol_pl)
#define MPP_pl (MPP_pl_amt_ / Vol_pl)
#define CLP_pl (CLP_pl_amt_ / Vol_pl)
#define MEP_pl (MEP_pl_amt_ / Vol_pl)
#define BFU_eE_pl (BFU_eE_pl_amt_ / Vol_pl)
#define BFU_lE_pl (BFU_lE_pl_amt_ / Vol_pl)
#define CMP_pl (CMP_pl_amt_ / Vol_pl)
#define GMP_pl (GMP_pl_amt_ / Vol_pl)
#define RET_pl (RET_pl_amt_ / Vol_pl)
#define ERY_pl (ERY_pl_amt_ / Vol_pl)
#define IL3_bm (IL3_bm_amt_ / Vol_bm)
#define IL3_pl (IL3_pl_amt_ / Vol_pl)
#define IL3_per (IL3_per_amt_ / Default)
#define SCF_bm (SCF_bm_amt_ / Vol_bm)
#define SCF_pl (SCF_pl_amt_ / Vol_pl)
#define scfR_pl (scfR_pl_amt_ / Vol_pl)
#define EPO_pl (EPO_pl_amt_ / Vol_pl)
#define EPO_bm (EPO_bm_amt_ / Vol_bm)
#define EPO_per (EPO_per_amt_ / Default)



$MAIN
if(TIME == 30000.0){	
  Dose_EPO_SC_IUkg = 300.0;
  EPO_pl_base = EPO_pl;
}
//if(TIME == 30000.0){	
//	EPO_pl_base = EPO_pl;
//  switch_inj = 1;
//}
//if(TIME == 30000.01){	
//  switch_inj = 0;
//}
//
switch_inj_0 = (0.0);
EPO_sc_0 = (0.0);
EPO_lym_0 = (0.0);
EPO_mRNA_amt__0 = (1.0) * Default;
HIF1a_0 = (351.0);
Roxa_pl_0 = (0.0);
Roxa_gut_0 = (0.0);
Roxa_per_0 = (0.0);
HSC_bm_amt__0 = (1428000.0) * Vol_bm;
MPP_bm_amt__0 = (1714000.0) * Vol_bm;
CMP_bm_amt__0 = (14283000.0) * Vol_bm;
CLP_bm_amt__0 = (1000.0) * Vol_bm;
MEP_bm_amt__0 = (35954188.0) * Vol_bm;
BFU_eE_bm_amt__0 = (12474255.0) * Vol_bm;
BFU_lE_bm_amt__0 = (16526076.0) * Vol_bm;
CFU_E_bm_amt__0 = (1000000.0) * Vol_bm;
PreproE_bm_amt__0 = (10000000.0) * Vol_bm;
ProE_bm_amt__0 = (20000000.0) * Vol_bm;
Baso_eE_bm_amt__0 = (41000000.0) * Vol_bm;
Baso_lE_bm_amt__0 = (86000000.0) * Vol_bm;
PolyE_bm_amt__0 = (158000000.0) * Vol_bm;
OrthoE_bm_amt__0 = (260000000.0) * Vol_bm;
RET_bm_amt__0 = (565000000.0) * Vol_bm;
GMP_bm_amt__0 = (17854000.0) * Vol_bm;
MKC_bm_amt__0 = (0.0) * Vol_bm;
NE_myelobl_bm_amt__0 = (25709000.0) * Vol_bm;
HSC_pl_amt__0 = (1000.0) * Vol_pl;
MPP_pl_amt__0 = (1000.0) * Vol_pl;
CLP_pl_amt__0 = (1000.0) * Vol_pl;
MEP_pl_amt__0 = (1000.0) * Vol_pl;
BFU_eE_pl_amt__0 = (1000.0) * Vol_pl;
BFU_lE_pl_amt__0 = (1000.0) * Vol_pl;
CMP_pl_amt__0 = (1000.0) * Vol_pl;
GMP_pl_amt__0 = (1000.0) * Vol_pl;
RET_pl_amt__0 = (93577981.0) * Vol_pl;
ERY_pl_amt__0 = (9357798165.0) * Vol_pl;
IL3_bm_amt__0 = (0.0) * Vol_bm;
IL3_pl_amt__0 = (0.0) * Vol_pl;
IL3_per_amt__0 = (0.0) * Default;
SCF_bm_amt__0 = (0.0) * Vol_bm;
SCF_pl_amt__0 = (0.0) * Vol_pl;
scfR_pl_amt__0 = (0.0) * Vol_pl;
R_SCF_bm_0 = (0.0);
scfR_bm_0 = (8000000.0);
EPO_pl_amt__0 = (1.333) * Vol_pl;
EPO_bm_amt__0 = (1.333) * Vol_bm;
EPO_per_amt__0 = (0.0) * Default;
epoR_bfu_le_0 = (80.0);
epoR_cfu_e_0 = (1100.0);
epoR_preproe_0 = (1100.0);
epoR_proe_0 = (1100.0);
epoR_baso_ee_0 = (600.0);
epoR_baso_le_0 = (600.0);
epoR_polye_0 = (320.0);
epoR_orthoe_0 = (80.0);
R_EPO_bfu_le_0 = (0.0);
R_EPO_cfu_e_0 = (0.0);
R_EPO_preproe_0 = (0.0);
R_EPO_proe_0 = (0.0);
R_EPO_baso_ee_0 = (0.0);
R_EPO_baso_le_0 = (0.0);
R_EPO_polye_0 = (0.0);
R_EPO_orthoe_0 = (0.0);

$ODE
// @Record ''
double EPO_pl_IU = EPO_pl * IU_per_pM;
// @Record ''
double EPO_pl_IU_corr = (EPO_pl - EPO_pl_base) * IU_per_pM;
// @Record ''
double EPO_bm_IU = EPO_bm * IU_per_pM;
// @Record ''
double EPO_cen_IU = EPO_pl * IU_per_pM;
// @Reaction ''
double V_syn_hif1a = Default * k_syn_hif1a;
// @Record ''
double Vmax_deg_hif1a = k_syn_hif1a / HIF1a_base;
// @Record ''
double Vol_bl = ((0.3669 * pow(BH, 3.0)) + (0.03219 * BW) + 0.6041) * sex + ((0.3561 * pow(BH, 3.0)) + (0.03308 * BW) + 0.1833) * (1.0 - sex);
// @Record ''
double HCT = 0.455 * sex + 0.402 * (1.0 - sex);
// @Compartment ''
double Vol_pl = Vol_bl * (1.0 - HCT);
// @Record ''
double Roxa_pl_uM = Roxa_pl / Mr_Roxa / Vol_pl * 1000.0;
// @Record ''
double HB = (RET_pl + ERY_pl) / (RET_pl_base + RBC_pl_base) * HB_base;
// @Reaction ''
double V_deg_hif1a = Default * Vmax_deg_hif1a * HIF1a * (1.0 - Roxa_pl_uM / (IC50_roxa + Roxa_pl_uM)) * pow((HB / HB_base), gamma_hb);
// @Reaction ''
double V_syn_epo_mrna = Default * k_syn_epo_mrna * pow((HIF1a / HIF1a_base), gamma_hif1a);
// @Record ''
double k_deg_epo_mrna = k_syn_epo_mrna;
// @Reaction ''
double V_deg_epo_mrna = Default * k_deg_epo_mrna * EPO_mRNA;
// @Reaction ''
double V_distr_epo_sc_lym = Default * k_epo_sc_lym * EPO_sc;
// @Reaction ''
double V_distr_epo_sc_pl = Default * k_epo_sc_pl * EPO_sc;
// @Reaction ''
double V_distr_epo_lym_pl = Default * k_epo_lym_pl * EPO_lym;
// @Reaction ''
double V_distr_epo_pl_per = Vol_pl * k_epo_pl_per * EPO_pl - Default * k_epo_per_pl * EPO_per;
// @Reaction ''
double V_distr_epo_pl_bm = Vol_pl * k_epo_pl_bm * EPO_pl - Vol_bm * k_epo_bm_pl * EPO_bm;
// @Reaction ''
double V_syn_epo_pl = Vol_pl * kbase_syn_epo_pl * EPO_mRNA;
// @Reaction ''
double V_deg_epo_pl = Vol_pl * (k_deg_epo_pl * EPO_pl + Vmax_deg_epo_pl * EPO_pl / (EPO_pl + Km_deg_epo_pl));
// @Reaction ''
double V_bind_epoR_bfu_le = k_dis_epoR * (EPO_bm * epoR_bfu_le / Kd_epoR - R_EPO_bfu_le);
// @Reaction ''
double V_bind_epoR_cfu_e = k_dis_epoR * (EPO_bm * epoR_cfu_e / Kd_epoR - R_EPO_cfu_e);
// @Reaction ''
double V_bind_epoR_preproe = k_dis_epoR * (EPO_bm * epoR_preproe / Kd_epoR - R_EPO_preproe);
// @Reaction ''
double V_bind_epoR_proe = k_dis_epoR * (EPO_bm * epoR_proe / Kd_epoR - R_EPO_proe);
// @Reaction ''
double V_bind_epoR_baso_ee = k_dis_epoR * (EPO_bm * epoR_baso_ee / Kd_epoR - R_EPO_baso_ee);
// @Reaction ''
double V_bind_epoR_baso_le = k_dis_epoR * (EPO_bm * epoR_baso_le / Kd_epoR - R_EPO_baso_le);
// @Reaction ''
double V_bind_epoR_polye = k_dis_epoR * (EPO_bm * epoR_polye / Kd_epoR - R_EPO_polye);
// @Reaction ''
double V_bind_epoR_orthoe = k_dis_epoR * (EPO_bm * epoR_orthoe / Kd_epoR - R_EPO_orthoe);
// @Reaction ''
double V_bind_epo_bm = Vol_bm / N_A_pmole * (V_bind_epoR_bfu_le * BFU_lE_bm + V_bind_epoR_cfu_e * CFU_E_bm + V_bind_epoR_preproe * PreproE_bm + V_bind_epoR_proe * ProE_bm + V_bind_epoR_baso_ee * Baso_eE_bm + V_bind_epoR_baso_le * Baso_lE_bm + V_bind_epoR_polye * PolyE_bm + V_bind_epoR_orthoe * OrthoE_bm) * 1000.0;
// @Record ''
double k_syn_epoR = k_syn_epoR_cfu_e / epoR_cfu_e_total_exp;
// @Reaction ''
double V_syn_epoR_bfu_le = k_syn_epoR * epoR_bfu_le_total_exp;
// @Reaction ''
double V_syn_epoR_cfu_e = k_syn_epoR * epoR_cfu_e_total_exp;
// @Reaction ''
double V_syn_epoR_preproe = k_syn_epoR * epoR_preproe_total_exp;
// @Reaction ''
double V_syn_epoR_proe = k_syn_epoR * epoR_proe_total_exp;
// @Reaction ''
double V_syn_epoR_baso_ee = k_syn_epoR * epoR_baso_ee_total_exp;
// @Reaction ''
double V_syn_epoR_baso_le = k_syn_epoR * epoR_baso_le_total_exp;
// @Reaction ''
double V_syn_epoR_polye = k_syn_epoR * epoR_polye_total_exp;
// @Reaction ''
double V_syn_epoR_orthoe = k_syn_epoR * epoR_orthoe_total_exp;
// @Reaction ''
double V_deg_epoR_bfu_le = k_deg_epoR * epoR_bfu_le;
// @Reaction ''
double V_deg_epoR_cfu_e = k_deg_epoR * epoR_cfu_e;
// @Reaction ''
double V_deg_epoR_preproe = k_deg_epoR * epoR_preproe;
// @Reaction ''
double V_deg_epoR_proe = k_deg_epoR * epoR_proe;
// @Reaction ''
double V_deg_epoR_baso_ee = k_deg_epoR * epoR_baso_ee;
// @Reaction ''
double V_deg_epoR_baso_le = k_deg_epoR * epoR_baso_le;
// @Reaction ''
double V_deg_epoR_polye = k_deg_epoR * epoR_polye;
// @Reaction ''
double V_deg_epoR_orthoe = k_deg_epoR * epoR_orthoe;
// @Reaction ''
double V_int_epoR_bfu_le = k_int_epoR * R_EPO_bfu_le;
// @Reaction ''
double V_int_epoR_cfu_e = k_int_epoR * R_EPO_cfu_e;
// @Reaction ''
double V_int_epoR_preproe = k_int_epoR * R_EPO_preproe;
// @Reaction ''
double V_int_epoR_proe = k_int_epoR * R_EPO_proe;
// @Reaction ''
double V_int_epoR_baso_ee = k_int_epoR * R_EPO_baso_ee;
// @Reaction ''
double V_int_epoR_baso_le = k_int_epoR * R_EPO_baso_le;
// @Reaction ''
double V_int_epoR_polye = k_int_epoR * R_EPO_polye;
// @Reaction ''
double V_int_epoR_orthoe = k_int_epoR * R_EPO_orthoe;
// @Reaction ''
double V_inj_epo_sc = switch_inj * Dose_EPO_SC_IUkg * BW / IU_per_pM / time_inj * F_EPO_SC;
// @Reaction ''
double V_inj_epo_iv = switch_inj * Dose_EPO_IV_IUkg * BW / IU_per_pM / time_inj;
// @Record ''
double AUC_EPO_pl_IU_h_L = AUC_EPO_pl / IU_per_pM / Vol_bm;
// @Record ''
double Dose_Roxa_PO = Dose_Roxa_PO_mg * F_Roxa;
// @Record ''
double Roxa_ug_L = Roxa_pl / Vol_pl * 1000.0;
// @Reaction ''
double V_absor_roxa_gut_pl = k_absorb_roxa * Roxa_gut;
// @Reaction ''
double V_distr_roxa_pl_per = k_roxa_pl_per * Roxa_pl - k_roxa_per_pl * Roxa_per;
// @Reaction ''
double V_deg_roxa_pl = k_deg_roxa_pl * Roxa_pl;
// @Record ''
double epoR_cfu_e_total = epoR_cfu_e + R_EPO_cfu_e;
// @Record ''
double epoR_occup = 100.0 * R_EPO_cfu_e / epoR_cfu_e_total;
// @Record ''
double portion_donation = (Vol_bl - donation) / Vol_bl;
// @Record ''
double cellularity = HSC_bm + MPP_bm + CMP_bm + CLP_bm + MEP_bm + BFU_eE_bm + BFU_lE_bm + CFU_E_bm + PreproE_bm + ProE_bm + Baso_eE_bm + Baso_lE_bm + PolyE_bm + OrthoE_bm + RET_bm + GMP_bm + MKC_bm + NE_myelobl_bm;
// @Record ''
double scfR_total = scfR_bm + R_SCF_bm;
// @Record ''
double coef_apo_scf = (scfR_bm + Imax_apo_scf * R_SCF_bm) / scfR_total;
// @Record ''
double coef_apo_il3 = 1.0 - Imax_apo_il3 * IL3_bm / IC50_apo_il3 / (1.0 + IL3_bm / IC50_apo_il3);
// @Record ''
double k_apo_hsc_bm = kbase_apo_hsc * coef_apo_scf * coef_apo_il3;
// @Reaction ''
double V_apo_hsc_bm = Vol_bm * k_apo_hsc_bm * HSC_bm;
// @Record ''
double k_apo_mpp_bm = kbase_apo_hsc * coef_apo_scf * coef_apo_il3;
// @Reaction ''
double V_apo_mpp_bm = Vol_bm * k_apo_mpp_bm * MPP_bm;
// @Reaction ''
double V_apo_clp_bm = Vol_bm * k_apo_clp_bm * CLP_bm;
// @Record ''
double k_apo_cmp_bm = kbase_apo_hsc * coef_apo_scf * coef_apo_il3;
// @Reaction ''
double V_apo_cmp_bm = Vol_bm * k_apo_cmp_bm * CMP_bm;
// @Reaction ''
double V_apo_gmp_bm = Vol_bm * k_apo_gmp_bm * GMP_bm;
// @Record ''
double k_apo_mep_bm = kbase_apo_hsc * coef_apo_scf;
// @Reaction ''
double V_apo_mep_bm = Vol_bm * k_apo_mep_bm * MEP_bm;
// @Record ''
double k_apo_bfu_ee_bm = kbase_apo_hsc * coef_apo_scf;
// @Reaction ''
double V_apo_bfu_ee_bm = Vol_bm * k_apo_bfu_ee_bm * BFU_eE_bm;
// @Record ''
double epoR_bfu_le_total = epoR_bfu_le + R_EPO_bfu_le;
// @Record ''
double k_apo_bfu_le_bm = kbase_apo_hsc * (coef_apo_scf + ka_epo_bfu * (epoR_bfu_le + Imax_apo_epo * R_EPO_bfu_le) / epoR_bfu_le_total + k_scf_epo_apo_bfu * R_EPO_bfu_le * R_SCF_bm / epoR_bfu_le_total / scfR_total);
// @Reaction ''
double V_apo_bfu_le_bm = Vol_bm * k_apo_bfu_le_bm * BFU_lE_bm;
// @Record ''
double k_apo_cfu_e_bm = kbase_apo_cfu_e_bm * epoR_cfu_e / epoR_cfu_e_total / (epoR_cfu_e / epoR_cfu_e_total + IC50_apo_epo_cfu_e);
// @Reaction ''
double V_apo_cfu_e_bm = Vol_bm * k_apo_cfu_e_bm * CFU_E_bm;
// @Record ''
double epoR_preproe_total = epoR_preproe + R_EPO_preproe;
// @Record ''
double k_apo_preproe_bm = kbase_apo_bm * (epoR_preproe + Imax_apo_epo * R_EPO_preproe) / epoR_preproe_total;
// @Reaction ''
double V_apo_preproe_bm = Vol_bm * k_apo_preproe_bm * PreproE_bm;
// @Record ''
double epoR_proe_total = epoR_proe + R_EPO_proe;
// @Record ''
double k_apo_proe_bm = kbase_apo_bm * (epoR_proe + Imax_apo_epo * R_EPO_proe) / epoR_proe_total;
// @Reaction ''
double V_apo_proe_bm = Vol_bm * k_apo_proe_bm * ProE_bm;
// @Record ''
double epoR_baso_ee_total = epoR_baso_ee + R_EPO_baso_ee;
// @Record ''
double k_apo_baso_ee_bm = kbase_apo_bm * (epoR_baso_ee + Imax_apo_epo * R_EPO_baso_ee) / epoR_baso_ee_total;
// @Reaction ''
double V_apo_baso_ee_bm = Vol_bm * k_apo_baso_ee_bm * Baso_eE_bm;
// @Record ''
double epoR_baso_le_total = epoR_baso_le + R_EPO_baso_le;
// @Record ''
double k_apo_baso_le_bm = kbase_apo_bm * (epoR_baso_le + Imax_apo_epo * R_EPO_baso_le) / epoR_baso_le_total;
// @Reaction ''
double V_apo_baso_le_bm = Vol_bm * k_apo_baso_le_bm * Baso_lE_bm;
// @Record ''
double epoR_polye_total = epoR_polye + R_EPO_polye;
// @Record ''
double k_apo_polye_bm = kbase_apo_bm * (epoR_polye + Imax_apo_epo * R_EPO_polye) / epoR_polye_total;
// @Reaction ''
double V_apo_polye_bm = Vol_bm * k_apo_polye_bm * PolyE_bm;
// @Record ''
double epoR_orthoe_total = epoR_orthoe + R_EPO_orthoe;
// @Record ''
double k_apo_orthoe_bm = kbase_apo_bm * (epoR_orthoe + Imax_apo_epo * R_EPO_orthoe) / epoR_orthoe_total;
// @Reaction ''
double V_apo_orthoe_bm = Vol_bm * k_apo_orthoe_bm * OrthoE_bm;
// @Reaction ''
double V_apo_ret_bm = Vol_bm * k_apo_ret_bm * RET_bm;
// @Reaction ''
double V_apo_mkc_bm = Vol_bm * k_apo_mkc_bm * MKC_bm;
// @Reaction ''
double V_apo_ne_myelobl_bm = Vol_bm * k_apo_ne_myelobl_bm * NE_myelobl_bm;
// @Reaction ''
double V_apo_hsc_pl = Vol_pl * k_apo_pl * HSC_pl;
// @Reaction ''
double V_apo_mpp_pl = Vol_pl * k_apo_pl * MPP_pl;
// @Reaction ''
double V_apo_clp_pl = Vol_pl * k_apo_pl * CLP_pl;
// @Reaction ''
double V_apo_cmp_pl = Vol_pl * k_apo_pl * CMP_pl;
// @Reaction ''
double V_apo_gmp_pl = Vol_pl * k_apo_pl * GMP_pl;
// @Reaction ''
double V_apo_mep_pl = Vol_pl * k_apo_pl * MEP_pl;
// @Reaction ''
double V_apo_bfu_ee_pl = Vol_pl * k_apo_pl * BFU_eE_pl;
// @Reaction ''
double V_apo_bfu_le_pl = Vol_pl * k_apo_pl * BFU_lE_pl;
// @Reaction ''
double V_apo_ret_pl = Vol_pl * k_apo_ret_pl * RET_pl;
// @Reaction ''
double V_apo_ery_pl = Vol_pl * k_apo_ery_pl * ERY_pl;
// @Record ''
double coef_pro_scf = (scfR_bm + Emax_pro_scf * R_SCF_bm) / scfR_total;
// @Record ''
double coef_pro_il3 = 1.0 + Emax_pro_il3 * IL3_bm / EC50_pro_il3 / (1.0 + IL3_bm / EC50_pro_il3);
// @Record ''
double k_pro_hsc_bm = kbase_pro * coef_pro_scf * coef_pro_il3;
// @Reaction ''
double V_pro_hsc_bm = Vol_bm * k_pro_hsc_bm * HSC_bm * EC50_proe / (EC50_proe + HSC_bm);
// @Record ''
double k_pro_mpp_bm = kbase_pro * coef_pro_scf * coef_pro_il3;
// @Reaction ''
double V_pro_mpp_bm = Vol_bm * k_pro_mpp_bm * MPP_bm * EC50_proe / (EC50_proe + MPP_bm);
// @Reaction ''
double V_pro_clp_bm = Vol_bm * k_pro_clp_bm * CLP_bm * EC50_proe / (EC50_proe + CLP_bm);
// @Record ''
double k_pro_cmp_bm = kbase_pro * coef_pro_scf * coef_pro_il3;
// @Reaction ''
double V_pro_cmp_bm = Vol_bm * k_pro_cmp_bm * CMP_bm * EC50_proe / (EC50_proe + CMP_bm);
// @Record ''
double k_pro_gmp_bm = kbase_pro * coef_pro_scf;
// @Reaction ''
double V_pro_gmp_bm = Vol_bm * k_pro_gmp_bm * GMP_bm * EC50_proe / (EC50_proe + GMP_bm);
// @Record ''
double k_pro_mep_bm = kbase_pro * coef_pro_scf;
// @Reaction ''
double V_pro_mep_bm = Vol_bm * k_pro_mep_bm * MEP_bm * EC50_proe / (EC50_proe + MEP_bm);
// @Record ''
double k_pro_bfu_bm = kbase_pro_bfu * coef_pro_scf;
// @Reaction ''
double V_pro_bfu_ee_bm = Vol_bm * k_pro_bfu_bm * BFU_eE_bm * EC50_proe / (EC50_proe + BFU_eE_bm);
// @Reaction ''
double V_pro_bfu_le_bm = Vol_bm * k_pro_bfu_bm * BFU_lE_bm * EC50_proe / (EC50_proe + BFU_lE_bm);
// @Record ''
double k_pro_cfu_e_bm = kmax_pro_cfu_e_bm * R_EPO_cfu_e / epoR_cfu_e_total / (R_EPO_cfu_e / epoR_cfu_e_total + EC50_pro_epo_cfu_e);
// @Reaction ''
double V_pro_cfu_e_bm = Vol_bm * k_pro_cfu_e_bm * CFU_E_bm * EC50_proe / (EC50_proe + CFU_E_bm);
// @Reaction ''
double V_dif_hsc_mpp_bm = Vol_bm * k_dif_hsc_mpp_bm * HSC_bm;
// @Reaction ''
double V_dif_mpp_clp_bm = Vol_bm * k_dif_mpp_clp_bm * MPP_bm;
// @Reaction ''
double V_dif_mpp_cmp_bm = Vol_bm * k_dif_mpp_cmp_bm * MPP_bm;
// @Reaction ''
double V_dif_cmp_gmp_bm = Vol_bm * k_dif_cmp_gmp_bm * CMP_bm;
// @Reaction ''
double V_dif_cmp_mep_bm = Vol_bm * k_dif_cmp_mep_bm * CMP_bm;
// @Reaction ''
double V_dif_gmp_myelobl_bm = Vol_bm * k_dif_gmp_myelobl_bm * GMP_bm;
// @Reaction ''
double V_dif_mep_bfu_ee_bm = Vol_bm * k_dif_mep_bfu_ee_bm * MEP_bm;
// @Reaction ''
double V_dif_mep_mkc_bm = Vol_bm * k_dif_mep_mkc_bm * MEP_bm;
// @Reaction ''
double V_dif_bfu_ee_bfu_le_bm = Vol_bm * k_dif_bfu_ee_bfu_le_bm * BFU_eE_bm;
// @Record ''
double Emax_dif_epo_bfu_le = Emax_dif_epo - (Emax_dif_epo - 1.0) * exp(-pow((epoR_bfu_le_total / coef_dif_epo), h_dif_epo));
// @Record ''
double k_dif_bfu_le_cfu_e_bm = kbase_dif_bfu_le_cfu_e_bm * (epoR_bfu_le + Emax_dif_epo_bfu_le * R_EPO_bfu_le) / epoR_bfu_le_total;
// @Reaction ''
double V_dif_bfu_le_cfu_e_bm = Vol_bm * k_dif_bfu_le_cfu_e_bm * BFU_lE_bm;
// @Record ''
double k_dif_cfu_e_preproe_bm = kmax_dif_cfu_e_preproe_bm * R_EPO_cfu_e / epoR_cfu_e_total / (R_EPO_cfu_e / epoR_cfu_e_total + EC50_dif_epo_cfu_e);
// @Reaction ''
double V_dif_cfu_e_preproe_bm = Vol_bm * k_dif_cfu_e_preproe_bm * CFU_E_bm;
// @Record ''
double Emax_dif_epo_preproe = Emax_dif_epo - (Emax_dif_epo - 1.0) * exp(-pow((epoR_preproe_total / coef_dif_epo), h_dif_epo));
// @Record ''
double k_dif_preproe_proe_bm = kbase_dif_bm * (epoR_preproe + Emax_dif_epo_preproe * R_EPO_preproe) / epoR_preproe_total;
// @Reaction ''
double V_dif_preproe_proe_bm = Vol_bm * k_dif_preproe_proe_bm * PreproE_bm;
// @Record ''
double Emax_dif_epo_proe = Emax_dif_epo - (Emax_dif_epo - 1.0) * exp(-pow((epoR_proe_total / coef_dif_epo), h_dif_epo));
// @Record ''
double k_dif_proe_baso_ee_bm = kbase_dif_bm * (epoR_proe + Emax_dif_epo_proe * R_EPO_proe) / epoR_proe_total;
// @Reaction ''
double V_dif_proe_baso_ee_bm = Vol_bm * k_dif_proe_baso_ee_bm * ProE_bm;
// @Record ''
double Emax_dif_epo_baso_ee = Emax_dif_epo - (Emax_dif_epo - 1.0) * exp(-pow((epoR_baso_ee_total / coef_dif_epo), h_dif_epo));
// @Record ''
double k_dif_baso_ee_baso_le_bm = kbase_dif_bm * (epoR_baso_ee + Emax_dif_epo_baso_ee * R_EPO_baso_ee) / epoR_baso_ee_total;
// @Reaction ''
double V_dif_baso_ee_baso_le_bm = Vol_bm * k_dif_baso_ee_baso_le_bm * Baso_eE_bm;
// @Record ''
double Emax_dif_epo_baso_le = Emax_dif_epo - (Emax_dif_epo - 1.0) * exp(-pow((epoR_baso_le_total / coef_dif_epo), h_dif_epo));
// @Record ''
double k_dif_baso_le_polye_bm = kbase_dif_bm * (epoR_baso_le + Emax_dif_epo_baso_le * R_EPO_baso_le) / epoR_baso_le_total;
// @Reaction ''
double V_dif_baso_le_polye_bm = Vol_bm * k_dif_baso_le_polye_bm * Baso_lE_bm;
// @Record ''
double Emax_dif_epo_polye = Emax_dif_epo - (Emax_dif_epo - 1.0) * exp(-pow((epoR_polye_total / coef_dif_epo), h_dif_epo));
// @Record ''
double k_dif_polye_orthoe_bm = kbase_dif_bm * (epoR_polye + Emax_dif_epo_polye * R_EPO_polye) / epoR_polye_total;
// @Reaction ''
double V_dif_polye_orthoe_bm = Vol_bm * k_dif_polye_orthoe_bm * PolyE_bm;
// @Record ''
double Emax_dif_epo_orthoe = Emax_dif_epo - (Emax_dif_epo - 1.0) * exp(-pow((epoR_orthoe_total / coef_dif_epo), h_dif_epo));
// @Record ''
double k_dif_orthoe_ret_bm = kbase_dif_orthoe_ret_bm * (epoR_orthoe + Emax_dif_epo_orthoe * R_EPO_orthoe) / epoR_orthoe_total;
// @Reaction ''
double V_dif_orthoe_ret_bm = Vol_bm * k_dif_orthoe_ret_bm * OrthoE_bm;
// @Reaction ''
double V_dif_ret_ery_pl = Vol_pl * k_dif_ret_ery_pl * RET_pl;
// @Reaction ''
double V_migr_hsc_bm_pl = Vol_bm * k_migr_bm_pl * HSC_bm;
// @Reaction ''
double V_migr_mpp_bm_pl = Vol_bm * k_migr_bm_pl * MPP_bm;
// @Reaction ''
double V_migr_clp_bm_pl = Vol_bm * k_migr_bm_pl * CLP_bm;
// @Reaction ''
double V_migr_cmp_bm_pl = Vol_bm * k_migr_bm_pl * CMP_bm;
// @Reaction ''
double V_migr_gmp_bm_pl = Vol_bm * k_migr_bm_pl * GMP_bm;
// @Reaction ''
double V_migr_mep_bm_pl = Vol_bm * k_migr_bm_pl * MEP_bm;
// @Reaction ''
double V_migr_bfu_ee_bm_pl = Vol_bm * k_migr_bm_pl * BFU_eE_bm;
// @Reaction ''
double V_migr_bfu_le_bm_pl = Vol_bm * k_migr_bm_pl * BFU_lE_bm;
// @Reaction ''
double V_migr_ret_bm_pl = Vol_bm * k_migr_ret_bm_pl * RET_bm;
// @Reaction ''
double V_prod_scf_pl = Vol_pl * k_prod_scf_pl;
// @Reaction ''
double V_prod_scf_bm = Vol_bm * k_prod_scf_bm;
// @Reaction ''
double V_deg_scf_pl = Vol_pl * k_deg_scf_pl * SCF_pl;
// @Reaction ''
double V_deg_scf_bm = Vol_bm * k_deg_scf_bm * SCF_bm;
// @Reaction ''
double V_bind_scfR_bm = k_ass_scfR * SCF_bm * scfR_bm - k_dis_scfR * R_SCF_bm;
// @Reaction ''
double V_bind_scf_bm = Vol_bm / N_A_pmole * V_bind_scfR_bm * (BFU_lE_bm + BFU_eE_bm + MEP_bm + CMP_bm + MPP_bm + HSC_bm + GMP_bm);
// @Reaction ''
double V_distr_scf_pl_bm = Vol_pl * k_scf_pl_bm * SCF_pl - Vol_bm * k_scf_bm_pl * SCF_bm;
// @Reaction ''
double V_distr_scfR_pl_bm = PS_scfR_pl_bm * (scfR_pl - scfR_bm);
// @Reaction ''
double V_prod_scfR_pl = Vol_pl * k_prod_scfR_pl;
// @Reaction ''
double V_deg_scfR_pl = Vol_pl * k_deg_scfR_pl * scfR_pl;
// @Reaction ''
double V_syn_scfR_bm = k_syn_scfR_bm;
// @Reaction ''
double V_deg_scfR_bm = k_deg_scfR_bm * scfR_bm;
// @Reaction ''
double V_int_scfR_bm = k_int_scfR * R_SCF_bm;
// @Reaction ''
double V_prod_il3_pl = Vol_pl * k_prod_il3_pl;
// @Reaction ''
double V_prod_il3_bm = Vol_bm * k_prod_il3_bm;
// @Reaction ''
double V_deg_il3_pl = Vol_pl * IL3_pl * (k_deg_il3_pl + Vm_deg_il3_pl / (Km_deg_il3_pl + IL3_pl));
// @Reaction ''
double V_deg_il3_bm = Vol_bm * IL3_bm * kbase_deg_il3_bm * (1.0 + kcell_deg_il3_bm * (HSC_bm + MPP_bm + CMP_bm) / (HSC_bm_base + MPP_bm_base + CMP_bm_base));
// @Reaction ''
double V_distr_il3_pl_per = Vol_pl * k_il3_pl_per * IL3_pl - Default * k_il3_per_pl * IL3_per;
// @Reaction ''
double V_distr_il3_pl_bm = Vol_pl * k_il3_pl_bm * IL3_pl - Vol_bm * k_il3_bm_pl * IL3_bm;
// @Record ''
double SCF_dimer_pl = Ka_scf * SCF_pl * SCF_pl;
// @Record ''
double deltaHB_perc = (HB - HB_base) / HB_base * 100.0;
// @Record ''
double deltaEPO_perc = (EPO_pl - EPO_pl_base) / EPO_pl_base * 100.0;
// @Record ''
double deltaRBC_perc = (ERY_pl - RBC_pl_base) / RBC_pl_base * 100.0;
// @Record ''
double deltaRET_perc = (RET_pl - RET_pl_base) / RET_pl_base * 100.0;
// @Record ''
double RET_pl_norm = RET_pl / RET_pl_base;
// @Record ''
double BasoE_giga = Baso_eE_bm * 1e-6;
// @Record ''
double BasoL_giga = Baso_lE_bm * 1e-6;
// @Record ''
double BFU_E = BFU_eE_bm + BFU_lE_bm;
// @Record ''
double Ortho_giga = OrthoE_bm * 1e-6;
// @Record ''
double Poly_giga = PolyE_bm * 1e-6;
// @Record ''
double Prepro_giga = PreproE_bm * 1e-6;
// @Record ''
double Pro_giga = ProE_bm * 1e-6;
// @Record ''
double RBC_blood_giga = ERY_pl * (1.0 - HCT) * 1e-6;
// @Record ''
double RET_blood_giga = RET_pl * (1.0 - HCT) * 1e-6;
// @Record ''
double RET_bm_giga = RET_bm * 1e-6;
// @Record ''
double RET_pl_giga = RET_pl * 1e-6;
// @Record ''
double T_dif_bfu_le_cfu_e_bm = 1.0 / k_dif_bfu_le_cfu_e_bm;
// @Record ''
double T_dif_preproe_proe_bm = 1.0 / k_dif_preproe_proe_bm;
// @Record ''
double T_dif_proe_baso_ee_bm = 1.0 / k_dif_proe_baso_ee_bm;
// @Record ''
double T_dif_baso_ee_baso_le_bm = 1.0 / k_dif_baso_ee_baso_le_bm;
// @Record ''
double T_dif_baso_le_polye_bm = 1.0 / k_dif_baso_le_polye_bm;
// @Record ''
double T_dif_polye_orthoe_bm = 1.0 / k_dif_polye_orthoe_bm;
// @Record ''
double T_dif_orthoe_ret_bm = 1.0 / k_dif_orthoe_ret_bm;
// @Record ''
double T_dif_cfu_e_preproe_bm = 1.0 / (k_dif_cfu_e_preproe_bm + 0.01);

dxdt_switch_inj = 0;
dxdt_EPO_sc = (-1)*V_distr_epo_sc_lym + (-1)*V_distr_epo_sc_pl + (1)*V_inj_epo_sc;
dxdt_EPO_lym = (1)*V_distr_epo_sc_lym + (-1)*V_distr_epo_lym_pl;
dxdt_EPO_mRNA_amt_ = (1)*V_syn_epo_mrna + (-1)*V_deg_epo_mrna;
dxdt_HIF1a = (1)*V_syn_hif1a + (-1)*V_deg_hif1a;
dxdt_Roxa_pl = (1)*V_absor_roxa_gut_pl + (-1)*V_distr_roxa_pl_per + (-1)*V_deg_roxa_pl;
dxdt_Roxa_gut = (-1)*V_absor_roxa_gut_pl;
dxdt_Roxa_per = (1)*V_distr_roxa_pl_per;
dxdt_HSC_bm_amt_ = (-1)*V_apo_hsc_bm + (-1)*V_pro_hsc_bm + (2)*V_pro_hsc_bm + (-1)*V_dif_hsc_mpp_bm + (-1)*V_migr_hsc_bm_pl;
dxdt_MPP_bm_amt_ = (-1)*V_apo_mpp_bm + (-1)*V_pro_mpp_bm + (2)*V_pro_mpp_bm + (1)*V_dif_hsc_mpp_bm + (-1)*V_dif_mpp_clp_bm + (-1)*V_dif_mpp_cmp_bm + (-1)*V_migr_mpp_bm_pl;
dxdt_CMP_bm_amt_ = (-1)*V_apo_cmp_bm + (-1)*V_pro_cmp_bm + (2)*V_pro_cmp_bm + (1)*V_dif_mpp_cmp_bm + (-1)*V_dif_cmp_gmp_bm + (-1)*V_dif_cmp_mep_bm + (-1)*V_migr_cmp_bm_pl;
dxdt_CLP_bm_amt_ = (-1)*V_apo_clp_bm + (-1)*V_pro_clp_bm + (2)*V_pro_clp_bm + (1)*V_dif_mpp_clp_bm + (-1)*V_migr_clp_bm_pl;
dxdt_MEP_bm_amt_ = (-1)*V_apo_mep_bm + (-1)*V_pro_mep_bm + (2)*V_pro_mep_bm + (1)*V_dif_cmp_mep_bm + (-1)*V_dif_mep_bfu_ee_bm + (-1)*V_dif_mep_mkc_bm + (-1)*V_migr_mep_bm_pl;
dxdt_BFU_eE_bm_amt_ = (-1)*V_apo_bfu_ee_bm + (-1)*V_pro_bfu_ee_bm + (2)*V_pro_bfu_ee_bm + (1)*V_dif_mep_bfu_ee_bm + (-1)*V_dif_bfu_ee_bfu_le_bm + (-1)*V_migr_bfu_ee_bm_pl;
dxdt_BFU_lE_bm_amt_ = (-1)*V_apo_bfu_le_bm + (-1)*V_pro_bfu_le_bm + (2)*V_pro_bfu_le_bm + (2)*V_dif_bfu_ee_bfu_le_bm + (-1)*V_dif_bfu_le_cfu_e_bm + (-1)*V_migr_bfu_le_bm_pl;
dxdt_CFU_E_bm_amt_ = (-1)*V_apo_cfu_e_bm + (-1)*V_pro_cfu_e_bm + (2)*V_pro_cfu_e_bm + (2)*V_dif_bfu_le_cfu_e_bm + (-1)*V_dif_cfu_e_preproe_bm;
dxdt_PreproE_bm_amt_ = (-1)*V_apo_preproe_bm + (2)*V_dif_cfu_e_preproe_bm + (-1)*V_dif_preproe_proe_bm;
dxdt_ProE_bm_amt_ = (-1)*V_apo_proe_bm + (2)*V_dif_preproe_proe_bm + (-1)*V_dif_proe_baso_ee_bm;
dxdt_Baso_eE_bm_amt_ = (-1)*V_apo_baso_ee_bm + (2)*V_dif_proe_baso_ee_bm + (-1)*V_dif_baso_ee_baso_le_bm;
dxdt_Baso_lE_bm_amt_ = (-1)*V_apo_baso_le_bm + (2)*V_dif_baso_ee_baso_le_bm + (-1)*V_dif_baso_le_polye_bm;
dxdt_PolyE_bm_amt_ = (-1)*V_apo_polye_bm + (2)*V_dif_baso_le_polye_bm + (-1)*V_dif_polye_orthoe_bm;
dxdt_OrthoE_bm_amt_ = (-1)*V_apo_orthoe_bm + (2)*V_dif_polye_orthoe_bm + (-1)*V_dif_orthoe_ret_bm;
dxdt_RET_bm_amt_ = (-1)*V_apo_ret_bm + (1)*V_dif_orthoe_ret_bm + (-1)*V_migr_ret_bm_pl;
dxdt_GMP_bm_amt_ = (-1)*V_apo_gmp_bm + (-1)*V_pro_gmp_bm + (2)*V_pro_gmp_bm + (1)*V_dif_cmp_gmp_bm + (-1)*V_dif_gmp_myelobl_bm + (-1)*V_migr_gmp_bm_pl;
dxdt_MKC_bm_amt_ = (-1)*V_apo_mkc_bm + (1)*V_dif_mep_mkc_bm;
dxdt_NE_myelobl_bm_amt_ = (-1)*V_apo_ne_myelobl_bm + (1)*V_dif_gmp_myelobl_bm;
dxdt_HSC_pl_amt_ = (-1)*V_apo_hsc_pl + (1)*V_migr_hsc_bm_pl;
dxdt_MPP_pl_amt_ = (-1)*V_apo_mpp_pl + (1)*V_migr_mpp_bm_pl;
dxdt_CLP_pl_amt_ = (-1)*V_apo_clp_pl + (1)*V_migr_clp_bm_pl;
dxdt_MEP_pl_amt_ = (-1)*V_apo_mep_pl + (1)*V_migr_mep_bm_pl;
dxdt_BFU_eE_pl_amt_ = (-1)*V_apo_bfu_ee_pl + (1)*V_migr_bfu_ee_bm_pl;
dxdt_BFU_lE_pl_amt_ = (-1)*V_apo_bfu_le_pl + (1)*V_migr_bfu_le_bm_pl;
dxdt_CMP_pl_amt_ = (-1)*V_apo_cmp_pl + (1)*V_migr_cmp_bm_pl;
dxdt_GMP_pl_amt_ = (-1)*V_apo_gmp_pl + (1)*V_migr_gmp_bm_pl;
dxdt_RET_pl_amt_ = (-1)*V_apo_ret_pl + (-1)*V_dif_ret_ery_pl + (1)*V_migr_ret_bm_pl;
dxdt_ERY_pl_amt_ = (-1)*V_apo_ery_pl + (1)*V_dif_ret_ery_pl;
dxdt_IL3_bm_amt_ = (1)*V_prod_il3_bm + (-1)*V_deg_il3_bm + (1)*V_distr_il3_pl_bm;
dxdt_IL3_pl_amt_ = (1)*V_prod_il3_pl + (-1)*V_deg_il3_pl + (-1)*V_distr_il3_pl_per + (-1)*V_distr_il3_pl_bm;
dxdt_IL3_per_amt_ = (1)*V_distr_il3_pl_per;
dxdt_SCF_bm_amt_ = (1)*V_prod_scf_bm + (-1)*V_deg_scf_bm + (-1)*V_bind_scf_bm + (1)*V_distr_scf_pl_bm;
dxdt_SCF_pl_amt_ = (1)*V_prod_scf_pl + (-1)*V_deg_scf_pl + (-1)*V_distr_scf_pl_bm;
dxdt_scfR_pl_amt_ = (-1)*V_distr_scfR_pl_bm + (1)*V_prod_scfR_pl + (-1)*V_deg_scfR_pl;
dxdt_R_SCF_bm = (-1)*V_int_scfR_bm + (1)*V_bind_scfR_bm;
dxdt_scfR_bm = (1)*V_distr_scfR_pl_bm + (1)*V_syn_scfR_bm + (-1)*V_deg_scfR_bm + (-1)*V_bind_scfR_bm;
dxdt_EPO_pl_amt_ = (1)*V_distr_epo_sc_pl + (1)*V_distr_epo_lym_pl + (-1)*V_distr_epo_pl_per + (-1)*V_distr_epo_pl_bm + (1)*V_syn_epo_pl + (-1)*V_deg_epo_pl + (1)*V_inj_epo_iv;
dxdt_EPO_bm_amt_ = (1)*V_distr_epo_pl_bm + (-1)*V_bind_epo_bm;
dxdt_EPO_per_amt_ = (1)*V_distr_epo_pl_per;
dxdt_epoR_bfu_le = (1)*V_syn_epoR_bfu_le + (-1)*V_deg_epoR_bfu_le + (-1)*V_bind_epoR_bfu_le;
dxdt_epoR_cfu_e = (1)*V_syn_epoR_cfu_e + (-1)*V_deg_epoR_cfu_e + (-1)*V_bind_epoR_cfu_e;
dxdt_epoR_preproe = (1)*V_syn_epoR_preproe + (-1)*V_deg_epoR_preproe + (-1)*V_bind_epoR_preproe;
dxdt_epoR_proe = (1)*V_syn_epoR_proe + (-1)*V_deg_epoR_proe + (-1)*V_bind_epoR_proe;
dxdt_epoR_baso_ee = (1)*V_syn_epoR_baso_ee + (-1)*V_deg_epoR_baso_ee + (-1)*V_bind_epoR_baso_ee;
dxdt_epoR_baso_le = (1)*V_syn_epoR_baso_le + (-1)*V_deg_epoR_baso_le + (-1)*V_bind_epoR_baso_le;
dxdt_epoR_polye = (1)*V_syn_epoR_polye + (-1)*V_deg_epoR_polye + (-1)*V_bind_epoR_polye;
dxdt_epoR_orthoe = (1)*V_syn_epoR_orthoe + (-1)*V_deg_epoR_orthoe + (-1)*V_bind_epoR_orthoe;
dxdt_R_EPO_bfu_le = (-1)*V_int_epoR_bfu_le + (1)*V_bind_epoR_bfu_le;
dxdt_R_EPO_cfu_e = (-1)*V_int_epoR_cfu_e + (1)*V_bind_epoR_cfu_e;
dxdt_R_EPO_preproe = (-1)*V_int_epoR_preproe + (1)*V_bind_epoR_preproe;
dxdt_R_EPO_proe = (-1)*V_int_epoR_proe + (1)*V_bind_epoR_proe;
dxdt_R_EPO_baso_ee = (-1)*V_int_epoR_baso_ee + (1)*V_bind_epoR_baso_ee;
dxdt_R_EPO_baso_le = (-1)*V_int_epoR_baso_le + (1)*V_bind_epoR_baso_le;
dxdt_R_EPO_polye = (-1)*V_int_epoR_polye + (1)*V_bind_epoR_polye;
dxdt_R_EPO_orthoe = (-1)*V_int_epoR_orthoe + (1)*V_bind_epoR_orthoe;

$CAPTURE @annotated
EPO_pl_IU :  (UL/L)
EPO_pl_base :  (pM)
EPO_pl_IU_corr :  (UL/L)
EPO_bm_IU :  (UL/L)
EPO_cen_IU :  (UL/L)
//switch_inj :  (UL)
k_deg_epo_mrna :  (1/h)
Vmax_deg_hif1a :  (1/h)
Dose_EPO_SC_IUkg :  (UL/kg)
F_EPO_SC :  
Dose_EPO_IV_IUkg :  (UL/kg)
EPO_mRNA :  (pM)
k_syn_epoR :  (1/h)
AUC_EPO_pl :  (pM*h)
AUC_EPO_pl_IU_h_L :  (UL*h/L)
Dose_Roxa_PO :  (mg)
AUC_Roxa_pl :  (mg*h)
Roxa_pl_uM :  (uM)
Roxa_ug_L :  (ug/L)
Vol_bl :  (L)
HCT :  (UL)
RET_pl_base :  (kcell/L)
RBC_pl_base :  (kcell/L)
HB_base :  (g/dL)
donation :  (L)
Vol_pl :  (L)
Vol_bm :  (L)
Vol_sp :  (L)
Default :  (L)
HSC_bm :  (kcell/L)
MPP_bm :  (kcell/L)
CMP_bm :  (kcell/L)
CLP_bm :  (kcell/L)
MEP_bm :  (kcell/L)
BFU_eE_bm :  (kcell/L)
BFU_lE_bm :  (kcell/L)
CFU_E_bm :  (kcell/L)
PreproE_bm :  (kcell/L)
ProE_bm :  (kcell/L)
Baso_eE_bm :  (kcell/L)
Baso_lE_bm :  (kcell/L)
PolyE_bm :  (kcell/L)
OrthoE_bm :  (kcell/L)
RET_bm :  (kcell/L)
GMP_bm :  (kcell/L)
MKC_bm :  (kcell/L)
NE_myelobl_bm :  (kcell/L)
HSC_pl :  (kcell/L)
MPP_pl :  (kcell/L)
CLP_pl :  (kcell/L)
MEP_pl :  (kcell/L)
BFU_eE_pl :  (kcell/L)
BFU_lE_pl :  (kcell/L)
CMP_pl :  (kcell/L)
GMP_pl :  (kcell/L)
RET_pl :  (kcell/L)
ERY_pl :  (kcell/L)
IL3_bm :  (pM)
IL3_pl :  (pM)
IL3_per :  (pM)
SCF_bm :  (pM)
SCF_pl :  (pM)
scfRs_bm :  (pM)
scfR_pl :  (pM)
EPO_pl :  (pM)
EPO_bm :  (pM)
EPO_per :  (pM)
epoR_baso_ee_total :  (item/cell)
epoR_baso_le_total :  (item/cell)
epoR_bfu_le_total :  (item/cell)
epoR_cfu_e_total :  (item/cell)
epoR_orthoe_total :  (item/cell)
epoR_polye_total :  (item/cell)
epoR_preproe_total :  (item/cell)
epoR_proe_total :  (item/cell)
epoR_occup :  (percent)
scfR_total :  (item/kcell)
coef_pro_scf :  (item/item)
coef_pro_il3 :  (UL)
k_pro_hsc_bm :  (1/h)
k_pro_mpp_bm :  (1/h)
k_pro_cmp_bm :  (1/h)
k_pro_gmp_bm :  (1/h)
k_pro_mep_bm :  (1/h)
k_pro_bfu_bm :  (1/h)
k_pro_cfu_e_bm :  (1/h)
Emax_dif_epo_baso_ee :  (UL)
Emax_dif_epo_baso_le :  (UL)
Emax_dif_epo_bfu_le :  (UL)
Emax_dif_epo_orthoe :  (UL)
Emax_dif_epo_polye :  (UL)
Emax_dif_epo_preproe :  (UL)
Emax_dif_epo_proe :  (UL)
k_dif_bfu_le_cfu_e_bm :  (1/h)
k_dif_preproe_proe_bm :  (1/h)
k_dif_proe_baso_ee_bm :  (1/h)
k_dif_baso_ee_baso_le_bm :  (1/h)
k_dif_baso_le_polye_bm :  (1/h)
k_dif_polye_orthoe_bm :  (1/h)
k_dif_orthoe_ret_bm :  (1/h)
k_dif_cfu_e_preproe_bm :  (1/h)
coef_apo_scf :  (item/item)
coef_apo_il3 :  (UL)
k_apo_hsc_bm :  (1/h)
k_apo_mpp_bm :  (1/h)
k_apo_cmp_bm :  (1/h)
k_apo_mep_bm :  (1/h)
k_apo_bfu_ee_bm :  (1/h)
k_apo_bfu_le_bm :  (1/h)
k_apo_cfu_e_bm :  (1/h)
k_apo_preproe_bm :  (1/h)
k_apo_proe_bm :  (1/h)
k_apo_baso_ee_bm :  (1/h)
k_apo_baso_le_bm :  (1/h)
k_apo_polye_bm :  (1/h)
k_apo_orthoe_bm :  (1/h)
portion_donation :  (L/L)
cellularity :  (kcell/L)
SCF_dimer_pl :  
HB :  
deltaHB_perc :  (percent)
deltaEPO_perc :  (percent)
deltaRBC_perc :  (percent)
deltaRET_perc :  (percent)
RET_pl_norm :  (kcell/kcell)
BasoE_giga :  
BasoL_giga :  
BFU_E :  
Ortho_giga :  
Poly_giga :  
Prepro_giga :  
Pro_giga :  
RBC_blood_giga :  
RET_blood_giga :  
RET_bm_giga :  
RET_pl_giga :  
T_dif_bfu_le_cfu_e_bm :  (h)
T_dif_preproe_proe_bm :  (h)
T_dif_proe_baso_ee_bm :  (h)
T_dif_baso_ee_baso_le_bm :  (h)
T_dif_baso_le_polye_bm :  (h)
T_dif_polye_orthoe_bm :  (h)
T_dif_orthoe_ret_bm :  (h)
T_dif_cfu_e_preproe_bm :  (h)
V_inj_epo_sc:
EPO_sc_0:
//
$PREAMBLE
double EPO_pl_base = 1.1;
//double switch_inj = 0.0;
double Dose_EPO_SC_IUkg = 0.0;
//double Dose_EPO_SC_IUkg = 300.0;
double F_EPO_SC = 1.0;
double Dose_EPO_IV_IUkg = 0.0;
//double EPO_sc = 0.0;
//double EPO_lym = 0.0;
double Default = 1.0;
//double EPO_mRNA = 1.0;
//double HIF1a = 351.0;
double AUC_EPO_pl = 0.0;
//double Roxa_pl = 0.0;
//double Roxa_gut = 0.0;
//double Roxa_per = 0.0;
double AUC_Roxa_pl = 0.0;
double RET_pl_base = 93577981.0;
double RBC_pl_base = 9357798165.0;
double HB_base = 13.5;
double donation = 0.0;
double Vol_sp = 0.5024;
double Vol_bm = 0.5024;
//double HSC_bm = 1428000.0;
//double MPP_bm = 1714000.0;
//double CMP_bm = 14283000.0;
//double CLP_bm = 1000.0;
//double MEP_bm = 35954188.0;
//double BFU_eE_bm = 12474255.0;
//double BFU_lE_bm = 16526076.0;
//double CFU_E_bm = 1000000.0;
//double PreproE_bm = 10000000.0;
//double ProE_bm = 20000000.0;
//double Baso_eE_bm = 41000000.0;
//double Baso_lE_bm = 86000000.0;
//double PolyE_bm = 158000000.0;
//double OrthoE_bm = 260000000.0;
//double RET_bm = 565000000.0;
//double GMP_bm = 17854000.0;
//double MKC_bm = 0.0;
//double NE_myelobl_bm = 25709000.0;
//double HSC_pl = 1000.0;
//double MPP_pl = 1000.0;
//double CLP_pl = 1000.0;
//double MEP_pl = 1000.0;
//double BFU_eE_pl = 1000.0;
//double BFU_lE_pl = 1000.0;
//double CMP_pl = 1000.0;
//double GMP_pl = 1000.0;
//double RET_pl = 93577981.0;
//double ERY_pl = 9357798165.0;
//double IL3_bm = 0.0;
//double IL3_pl = 0.0;
//double IL3_per = 0.0;
//double SCF_bm = 0.0;
//double SCF_pl = 0.0;
double scfRs_bm = 0.0;
//double scfR_pl = 0.0;
//double R_SCF_bm = 0.0;
//double scfR_bm = 8000000.0;
//double EPO_pl = 1.333;
//double EPO_bm = 1.333;
//double EPO_per = 0.0;
//double epoR_bfu_le = 80.0;
//double epoR_cfu_e = 1100.0;
//double epoR_preproe = 1100.0;
//double epoR_proe = 1100.0;
//double epoR_baso_ee = 600.0;
//double epoR_baso_le = 600.0;
//double epoR_polye = 320.0;
//double epoR_orthoe = 80.0;
//double R_EPO_bfu_le = 0.0;
//double R_EPO_cfu_e = 0.0;
//double R_EPO_preproe = 0.0;
//double R_EPO_proe = 0.0;
//double R_EPO_baso_ee = 0.0;
//double R_EPO_baso_le = 0.0;
//double R_EPO_polye = 0.0;
//double R_EPO_orthoe = 0.0;