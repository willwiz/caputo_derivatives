#!/usr/bin/env pathon3
import matplotlib.pyplot as plt
import numpy as np



def bar_plot(fout, data, label):
  label = np.asanyarray(label)
  # set width of bar
  barWidth = 0.33

  fig, axs = plt.subplots(2,2,figsize=(10,6), dpi=600, gridspec_kw={'width_ratios': [11, 2]})

  indices = [[i for i in range(1,11)], [i for i in range(11,21)], [i for i in range(21,23)]]

  pos = [(0,0), (1,0), (1,1)]
  for (i,j), k in zip(pos, indices):
    # Set position of bar on X axis
    br1 = np.arange(len(k)) + (0.5 - barWidth)
    br2 = [x + barWidth for x in br1]
    # Make the plot
    axs[i,j].bar(br1, 100* data[k, 0], color ='r', width = barWidth,
            edgecolor ='grey', label ='ZZ')
    axs[i,j].bar(br2, 100* data[k, 3], color ='b', width = barWidth,
            edgecolor ='grey', label ='CC')
    axs[i,j].set_xticks([r + barWidth for r in range(len(br1))])
    axs[i,j].set_xticklabels(label[k])
  # Adding Xticks
  axs[1,0].set_xlabel('Cycling', fontweight ='bold', fontsize = 15)
  axs[1,1].set_xlabel('Relaxation',  fontweight ='bold', fontsize = 15)
  axs[0,0].set_ylabel('Absolute err (%)', fontweight ='bold', fontsize = 15)
  axs[1,0].set_ylabel('Absolute err (%)', fontweight ='bold', fontsize = 15)
  fig.delaxes(axs[0,1])

  handles, labels = axs[0,0].get_legend_handles_labels()
  fig.legend(handles, labels, loc='upper right', prop={'size': 15},
    framealpha=1)

  plt.tight_layout()
  plt.savefig(fout, bbox_inches='tight', transparent=True)
  plt.close()