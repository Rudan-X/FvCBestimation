library(DiagrammeR)

grViz("
digraph {
  graph [layout = dot, rankdir = LR]

  node [shape=box style=filled fontname=Helvetica]

  ## Level 2: Models (gray)
  subgraph cluster_level2 {
    FvCB98 [label='Farquhar (1980)']
    Collatz92  [label='Collatz (1992)']
    Ethier04   [label='Ethier &amp; Livingston (2004)']
    Caemmerer00[label='von Caemmerer (2000)']
    Yin04      [label='Yin (2004)']
    Price11    [label='Price (2011)']
    Tholen12   [label='Tholen (2012)']
    Evans12    [label='Evans &amp; von Caemmerer (2012)']
    Busch18    [label='Busch (2018)']
    Busch20    [label='Busch (2020)']
  }

############# Estimation methods ####################

  ## Level 3: Estimation methods (purple-gray) ======
  subgraph cluster_level4 {
    Harley92  [label='Harley (1992): const & variable J' fillcolor='#E5D4FF']
    Dubois07  [label='Dubois (2007): Nonlinear Least Squares (SAS)' fillcolor='#D9E2EC']
    Sharkey07  [label='Sharkey (2007): NLS (Excel)' fillcolor='#FFF1B6']
    Su09  [label='Su (2009): Genetic algorithm (MATLAB)' fillcolor='#D9E2EC']
    Patrick09  [label='Patrick (2009): Hierarchical Bayesian (WinBUGS)' fillcolor='#D9E2EC']
    Yin09a  [label='Yin (2009): NLS (SAS)' fillcolor='#E5D4FF']
    Gu10  [label='Gu (2010): Hybrid of gradient and non-gradient approaches (LeafWeb)' fillcolor='#D9E2EC']
    Yin09b  [label='Yin (2011): NLS (SAS)' fillcolor='#E6FFFB']
    Moualeu16 [label='Moualeu (2016):  NLS (SAS)' fillcolor='#D9E2EC']
    Plantecophys  [label='Plantecophys (Duursma 2015): NLS' fillcolor='#D9E2EC']
    Xiao21  [label='Xiao (2021): Bayesian' fillcolor='#E6FFFB']
    msu21 [label='msuRACiFit (Gregory 2021): Excel' fillcolor='#FFF1B6']
    PhotoGEA [label='PhotoGEA (Lochocki 2025): NLS' fillcolor='#FFF1B6']
    PhoTorch [label='PhoTorch (Lei 2025): ADAM optimizer' fillcolor='#D9E2EC']

    Harley92_pars [label='gm (constant & variable J methods)'  fillcolor='#E5D4FF']
    Dubois07_pars [label='Vcmax, J, Rd'    fillcolor='#D9E2EC']
    Sharkey07_pars [label='Vcmax, J, Vtpu, Rd (gm)'    fillcolor='#FFF1B6']
    Su09_pars [label='Vcmax, J, Vtpu, Rd'    fillcolor='#D9E2EC']
    Patrick09_pars [label='Vcmax, Jmax, Rd, Kmc, Kmo, gm'    fillcolor='#D9E2EC']
    Yin09a_pars [label='gm (modified J, NRH-A, RH-A variants)'    fillcolor='#E5D4FF']
    Gu10_pars [label='Vcmax, J, Vtpu, Rd, alphaG, Kco, gamma_star, gm'    fillcolor='#D9E2EC']
    Yin09b_pars [label='Vcmax, Jmax, Vtpu, gm, Rd, Phi2(LL),s,k2(LL), theta'    fillcolor='#E6FFFB']
    Planteco_pars [label='Vcmax, Jmax,Rd'    fillcolor='#D9E2EC']
    Moualeu16_pars [label='Vcmax, J, Vtpu, Rd, s (gamma_star, Kmc, Kmo)'    fillcolor='#D9E2EC']
    Xiao21_pars [label='Vcmax, Jmax, Rd, gamma_star, s, Phi2(LL)'    fillcolor='#E6FFFB']
    msu_pars [label='Vcmax, J, Vtpu, Rd, gm, alphaG, alphaS'    fillcolor='#FFF1B6']
    Photogea_pars [label='Jmax, Vcmax, Vtpu, gm, Rd, alphaG, alphaS'    fillcolor='#FFF1B6']
    Photorch_pars [label='Jmax, Vcmax, Vtpu, alphaG'    fillcolor='#D9E2EC']
  }

  FvCB98 -> Harley92 [color='black' style=solid]
  FvCB98 -> Dubois07 [color='black' style=solid]
  Caemmerer00  -> Sharkey07 [color='black' style=solid]
  Ethier04  -> Sharkey07 [color='black' style=solid]
  Ethier04  -> Su09 [color='black' style=solid]
  Ethier04  -> Patrick09 [color='black' style=solid]
  Ethier04  -> Yin09a [color='black' style=solid]
  Caemmerer00  -> Gu10 [color='black' style=solid]
  Ethier04  -> Gu10 [color='black' style=solid]
  Yin04  -> Yin09b [color='black' style=solid]
  Caemmerer00  -> Xiao21 [color='black' style=solid]
  FvCB98 -> Moualeu16
  Busch18 -> msu21 [color='black' style=solid]
  Collatz92 -> Plantecophys [color='black' style=solid]
  Busch20 -> PhotoGEA [color='black' style=solid]
  Caemmerer00  -> PhoTorch [color='black' style=solid]
  # Moualeu16 -> PhotoGEA [color='black' style=solid]

  Harley92 -> Harley92_pars [label='A-Ci & CF (s3=2, r=4)',color='blue' style=dashed]
  Dubois07 -> Dubois07_pars [label='A-Ci (n=231, literature)',color='blue' style=dashed]
  Sharkey07 -> Sharkey07_pars [label='A-Ci',color='blue' style=dashed]
  Su09 -> Su09_pars [label='A-Ci (n=24, s3=21)',color='blue' style=dashed]
  Patrick09 -> Patrick09_pars [label='A–Ci (n=17) & A–Q (n=37), s3=4',color='blue' style=dashed]
  Yin09a -> Yin09a_pars [label='A-Ci & CF (at different N supply), r=4, literature',color='blue' style=dashed]
  Gu10 -> Gu10_pars [label='A-Ci (n=6, s3=3)',color='blue' style=dashed]
  Yin09b -> Yin09b_pars [label='A-Ci, A-Inc (up to 200) & CF (Cucumber)',color='blue' style=dashed]
  Moualeu16 -> Moualeu16_pars [label='',color='blue' style=dashed]
  Xiao21 -> Xiao21_pars [label='A-Ci (Inc=1500, 200, 100, 50)), CF',color='blue' style=dashed]
  Plantecophys -> Planteco_pars [label='A-Ci',color='blue' style=dashed]
  msu21 -> msu_pars [label='Data from Xiao 2021',color='blue' style=dashed]
  PhotoGEA -> Photogea_pars [label='GE + CF',color='blue' style=dashed]
  PhoTorch -> Photorch_pars [label='GE',color='blue' style=dashed]
}
")
