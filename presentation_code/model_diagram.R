library(DiagrammeR)

grViz("
digraph {
  graph [layout = dot, rankdir = LR]

  node [shape=box style=filled fontname=Helvetica]

  subgraph cluster_level1 { rank=same; style=invis;
    FvCB98_Ac   [label='FvCB98:Ac' fillcolor='#FDE2E4']
    FvCB98_Aj   [label='FvCB98:Aj' fillcolor='#E6FFFB']
    FvCB98 [label='Farquhar (1980)']
  }

  ## Level 2: Models (gray)
  subgraph cluster_level2 {
    Collatz92  [label='Collatz (1992)']
    Ethier04   [label='Ethier &amp; Livingston (2004)']
    Caemmerer00[label='von Caemmerer (2000)']
    Yin04      [label='Yin (2004)']
    Price11    [label='Price (2011)']
    Tholen12   [label='Tholen (2012)']
    Evans12    [label='Evans &amp; von Caemmerer (2012)']
    Busch18    [label='Busch (2018,2020)']

    Collatz92_An   [label='Collatz92:An'  fillcolor='white']
    Ethier04_Ac    [label='Ethier04:Ac'   fillcolor='#FDE2E4']
    Ethier04_Aj    [label='Ethier04:Aj'   fillcolor='#E6FFFB']
    Caemmerer00_Ap [label='Caemmerer00:Ap' fillcolor='#FFF1B6']
    Yin04_Aj       [label='Yin04:Aj'      fillcolor='#E6FFFB']
    Busch18_Aj     [label='Busch20:Aj'    fillcolor='#E6FFFB']
    Busch18_Ap     [label='Busch20:Ap'    fillcolor='#FFF1B6']
    Price11_gm     [label='Price11:gm' fillcolor='#E5D4FF']
    Evans12_gm     [label='New gm definition' fillcolor='#E5D4FF']
    Tholen12_gm    [label='Tholen12:gm' fillcolor='#E5D4FF']
  }


  # FvCB components to other models
  FvCB98_Ac -> FvCB98 [color='black' style=solid] #label='inherits'
  FvCB98_Aj -> FvCB98 [color='black' style=solid] #label='inherits'

  FvCB98_Ac -> Collatz92 [color='black' style=solid] #label='inherits'
  FvCB98_Aj -> Collatz92 [color='black' style=solid]

  FvCB98_Ac -> Ethier04 [color='black' style=solid]
  FvCB98_Aj -> Ethier04 [color='black' style=solid]

  FvCB98_Ac -> Yin04 [color='black' style=solid]

  FvCB98_Ac -> Busch18 [color='black' style=solid]

  # FvCV Model to model
  FvCB98 -> Caemmerer00 [color='black' style=solid]
  FvCB98 -> Price11 [color='black' style=solid]
  FvCB98 -> Tholen12 [color='black' style=solid]
  FvCB98 -> Evans12 [color='black' style=solid]


  # New model to new component
  Collatz92 -> Collatz92_An [label='Quadratic smoothing of An' color='blue' style=dashed]

  Caemmerer00 -> Caemmerer00_Ap [label='TPU from Sharkey' color='blue' style=dashed]

  Ethier04 -> Ethier04_Ac [label='Quadratic form' color='blue' style=dashed]
  Ethier04 -> Ethier04_Aj [label='Quadratic form' color='blue' style=dashed]

  Yin04 -> Yin04_Aj [label='Alternative ETC' color='blue' style=dashed]

  Price11 -> Price11_gm [label='gm & leakiness' color='blue' style=dashed]

  Tholen12 -> Tholen12_gm [label='New gm definition' color='blue' style=dashed]

  Evans12 -> Evans12_gm [label='gm & leakiness' color='blue' style=dashed]

  Busch18 -> Busch18_Aj [label='N assimilation' color='blue' style=dashed]
  Busch18 -> Busch18_Ap [color='blue' style=dashed]


}
")


ggsave("results/Figures/Model_diagram.png",width = 12, height = 7)
