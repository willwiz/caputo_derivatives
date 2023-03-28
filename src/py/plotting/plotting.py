#!/usr/bin/env python3
from cycler import cycler
import matplotlib.pyplot as plt
from itertools import cycle

plt.rc('lines', linewidth=0.75)
# plt.rc('axes', prop_cycle=(cycler('color', ['k', 'r', 'b', 'g', 'm', 'k']) +
#                            cycler('mec', ['k', 'r', 'b', 'g', 'm', 'k'])))
plt.rc('axes', prop_cycle=(cycler('color', ['k', 'r', 'b', 'g', 'm', 'k']) +
                           cycler('linestyle', ['None', '-', '-', '-', '-', '-'])+
                           cycler('marker', ['o','None', 'None', 'None', 'None', 'None'])+
                           cycler('mec', ['k', 'r', 'b', 'g', 'm', 'k'])))

def plot(fout, x, *ys, labelx="strain", labely=r'$S$ (kPa)', skip=1, msize=3,
    ylim=None, dpi=300, figsize=(4,3)):
  style = {'markersize': msize, 'markevery': skip, 'fillstyle': 'none', 'markeredgewidth': .2}
  fig, axs = plt.subplots(1, 1, dpi=dpi, figsize=figsize)
  axs.set(xlabel=labelx, ylabel=labely)
  for y in ys:
    axs.plot(x, y, **style)
  if ylim is not None:
    axs.set_ylim(ylim)
  fig.tight_layout()
  plt.savefig(fout, bbox_inches='tight', transparent=True)
  plt.close()

def plot_columns_vec(file_name, time, *datas, label='S', legend=['Raw','New'], ylim=None, dpi=300,
                components = {0:0, 1:3, 2:1, 3:2},
                linestyle = None,
                marker=None,
                style = {'markersize': 3, 'markevery': 1, 'fillstyle': 'none', 'markeredgewidth': .2}):
  if linestyle is None:
    linestyle = ['None', '-', '-', '-', '-', '-']
  if marker is None:
    marker = ['o','None', 'None', 'None', 'None', 'None']
  linecycler = cycle(linestyle)
  markercycler = cycle(marker)
  fig, axs = plt.subplots(len(components), 1, dpi=dpi, figsize=(8, 2*len(components)))
  for k, i in components.items():
    for m, d in enumerate(datas):
      axs[k].plot(time, d[:,i], label=legend[m], linestyle=next(linecycler), marker=next(markercycler), **style)
    # axs[k].set_title(fr'${label}_{{{i}}}$')
    if ylim is not None:
      axs[k].set_ylim(ylim)
  handles, labels = axs[0].get_legend_handles_labels()
  # fig.legend(handles, labels, loc='lower right',prop={'size': 11})
  fig.tight_layout()
  plt.savefig(file_name, bbox_inches='tight', transparent=True)
  plt.close()

def plot_rows_vec(file_name, time, *datas, label='S', legend=['Raw','New'], skip=1,
                msize=3, ylim=None, dpi=300,
                components = {0:0, 1:3, 2:1, 3:2}):
  style = {'markersize': msize, 'markevery': skip, 'fillstyle': 'none', 'markeredgewidth': .2}
  fig, axs = plt.subplots(1, len(components), dpi=dpi, figsize=(4*len(components), 3))
  for k, i in components.items():
    for m, d in enumerate(datas):
      axs[k].plot(time[:,i], d[:,i], label=legend[m], **style)
    # axs[k].set_title(fr'${label}_{{{i}}}$')
    if ylim is not None:
      axs[k].set_ylim(ylim)
  handles, labels = axs[0].get_legend_handles_labels()
  # fig.legend(handles, labels, loc='lower right',prop={'size': 11})
  fig.tight_layout()
  plt.savefig(file_name, bbox_inches='tight', transparent=True)
  plt.close()

def plot_tensor_time(file_name, time, *datas, label='S', legend=['Raw','New'], skip=1, msize=3,
    ylim=None, dpi=300, figsize=(8, 8), components = {0:(0,0), 1:(1,1), 2:(0,1), 3:(1,0)}):
  style = {'markersize': msize, 'markevery': skip, 'fillstyle': 'none', 'markeredgewidth': .2}
  fig, axs = plt.subplots(4, 1, dpi=dpi, figsize=figsize)
  for k, (i,j) in components.items():
    for m, d in enumerate(datas):
      axs[k].plot(time, d[:,i,j], label=legend[m], **style)
    # axs[k].set_title(fr'${label}_{{{i}{j}}}$')
    if ylim is not None:
      axs[k].set_ylim(ylim)
  handles, labels = axs[0].get_legend_handles_labels()
  # fig.legend(handles, labels, loc='lower right',prop={'size': 11})
  fig.tight_layout()
  plt.savefig(file_name, bbox_inches='tight', transparent=True)
  plt.close()