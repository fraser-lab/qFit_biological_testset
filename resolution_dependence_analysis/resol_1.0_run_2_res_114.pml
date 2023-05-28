# PyMOL script

# -------------------------
# General setup
# -------------------------

bg_color white;
space cmyk;
set orthoscopic, on;                # (default=off)
set valence, off;                   # (default=on)
set cartoon_side_chain_helper, on;
set cartoon_fancy_helices, on;
#set cartoon_highlight_color, grey90;
#set antialias, 2;

# Cartoon, w/ shading:
# (If the sum of these = 1, then perfect exposure. If > 1, then we lose brightness dynamic range (and the top end is full saturation).
set ray_trace_mode, 0;  ### THIS MUST BE RAY MODE 0 IF WE WANT TO USE RAY_VOLUME
set ray_shadow, off;  # Hard shadows, I don't like 'em!  # (default=on);
# set ray_volume, 1;  # Turn off ray volume if you want a transparent background!
set light_count, 8;
set ambient, 0.3;
set reflect, 0.4;
set direct, 0.8;  # 1.0 for very overexposed;
set specular, 0;
set ambient_occlusion_mode, 1;  #
set ambient_occlusion_smooth, 10;       # (default=10);
set ambient_occlusion_scale, 15;        # (default=25);

# DK in-house cartoon style
set cartoon_rect_length, 1.0;  # (default=1.40);  # (DKdefault=0.80);
set cartoon_oval_length, 1.0;  # (default=1.35);  # (DKdefault=0.80);
set stick_radius, 0.2;  # (default=0.25);
set solvent_radius, 1.6;  # smooths out single-water pockets in surface (default=1.4);
set sphere_scale, 0.15;
set dash_gap, 0.25;
set dash_color, black;
set mesh_width, 0.5;  # (default=1.0);
# set mesh_quality, 3;
# set mesh_grid_max, 1;
cartoon loop;
viewport 2000, 2000;

cmd.volume_ramp_new('density_vol_ramp', [\
    # level,     r,     g,     b, alph, \
      # L*a*b (100, -128, -128), \
      0.000, 0.000, 1.000, 1.000, 0.00, \
      # L*a*b ( 50, -128, -128), \
      1.000, 0.000, 0.663, 1.000, 0.10, \
      # L*a*b (  0, -128, -128), \
      2.000, 0.000, 0.252, 0.761, 0.05, \
      # L*a*b (  0, -128, -128), \
      2.000, 0.000, 0.252, 0.761, 0.00, \
]);

cmd.volume_ramp_new('difference_density_vol_ramp', [\
    # level,     r,     g,     b, alph, \
     -9.999, 1.000, 0.000, 0.000, 0.05, \
     -3.307, 1.000, 0.000, 0.000, 0.10, \
     -2.693, 1.000, 0.000, 0.000, 0.00, \
      2.693, 0.000, 1.000, 0.000, 0.00, \
      3.307, 0.000, 1.000, 0.000, 0.10, \
      9.999, 0.000, 1.000, 0.000, 0.05, \
]);

# -------------------------
# Fetch, create & organise object files
# -------------------------


resol="1.0"
run="2"
resid="114"
cmd.load("~/qfit_summit/summit_2/macdomain_maps/xray_noise_modeling/p1_box_pdbs/7kr0_shift_bfac_"+resol+"_p1.pdb")
cmd.load("~/qfit_summit/summit_2/macdomain_maps/resol_"+resol+"/run_"+run+"_final/qfit_bic/7kr0_resol_"+resol+"_shake_noisy_run_"+run+"_qFit.pdb")
cmd.load("/Users/ashraya/qfit_summit/summit_2/macdomain_maps/resol_"+resol+"/run_"+run+"_final/7kr0_resol_"+resol+"_shake_noisy_run_"+run+"_2mFo-DFc.ccp4")


sele m. 7kr0_resol_1.0_shake_noisy_run_2_qFit and i. 114;

isomesh resol_1.0.2fofc_mesh, 7kr0_resol_1.0_shake_noisy_run_2_2mFo-DFc, 1.0, selection=(sele), carve=2;

# -------------------------
# Set colours
# -------------------------
color 0x228833, m. 7kr0_shift_bfac_1.0_p1 and i. 114 and alt A;
color 0xCCBB44, m. 7kr0_shift_bfac_1.0_p1 and i. 114 and alt B;
color 0x66CCEE, m. 7kr0_resol_1.0_shake_noisy_run_2_qFit and i. 114 and alt B;
color 0xAA3377, m. 7kr0_resol_1.0_shake_noisy_run_2_qFit and i. 114 and alt C;

set mesh_color, density, *.2fofc_mesh
util.cnc;

# -------------------------
# Create view and render figures
# -------------------------

remove hydrogens;
hide everything;
show sticks, i. 114;
set stick_radius, 0.15;

set_view (\
     0.243568569,    0.121985838,   -0.961954832,\
     0.350380868,    0.914662600,    0.201364204,\
     0.905084968,   -0.385987282,    0.178387105,\
    -0.000000000,    0.000000000,  -22.215768814,\
    29.965373993,   15.873499870,   45.376750946,\
    19.380224228,   25.051309586,   20.000000000 )



show mesh, *fofc_mesh*;
orient sele

