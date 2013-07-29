%# file containing symmetry matrices:
syms =  load("Derek72grid/syms_para"); 

%# irreducible kgrid (wrt full group)
kirr = load("Derek72grid/k_irr"); 

%# assumes dirac bands are 4th and 5th col.
global energies = load("Derek72grid/en_shifted"); 

%# BGW kgrid:
global Nkx = 72; 

%# use wavefunctions and energies on shifted grid:
shift_flag = 1;
if(shift_flag == 1);
  %# unfolded shifted grid (for lda energies):
  global ks_shifted = load("Derek72grid/kshifted_full")(:,2:4);
  %# mapping from shifted grid to lda energies:
  global partner_shifted = load("Derek72grid/partner_shifted")
endif;

%# kpoint which is being refined:
global k = [0.03472222 0.68055556 0.00000000];

%# refined kgrid (new kgrid is Nkx*S(1) ):
S = [4 4]; 

%# freq. grid for BGW sigma calculation:
ws = [ -0.475246: 0.1: 0.524754 ]'; 

%# metric tensor in Gspace (atomic units):
global Gmat = [0.78013  -1.35122   0.00000;
               0.78013   1.35122   0.00000;
               0.00000   0.00000   0.33069];

%# load in output file from epsilon refinement
%# contains qcut and absolute Fermi energy
load("epsilon_output");

%# input for mini BZ plasmon pole model: 
%# file should contain mat, w_eV and int_eV 
%# can be obtained by running get_mBZ.m
load("Derek72grid/1eV_doping/dIepsI/dIepsI_dat");

%# list of qpoints that received epsilon-refinement:
qref_list = load("Derek72grid/1eV_doping/EpsInvDyn/qref_list"); 

%# directories for refined and unrefined heads of epsilon inverse:
ref_dir  = "Derek72grid/1eV_doping/EpsInvDyn/EpsInvDyn"; %# Nkx=72; 1 eV doping
unref_dir= "Derek72grid/1eV_doping/ParaUnrefined2/epsInvR"; %# Nkx=72; 1 eV doping

%# load symmetry reduced fine grid (has to be calculated first):
load_bzkf = 0; %# 0 means no

%# parameter for interpolation of heads of epsInv:
%# this freq. (in eV) should be between the carrier and the pi-plasmon peaks for all refined q
wc_eps = 1.975;

%# obtain the unrefined results (use with S=[1 1]):
global unref_flag = 0;
etaBGW = 0.15; %# broadening in unrefined epsilon calc