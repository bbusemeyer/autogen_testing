import scipy as sci
import mython as my
from random import gauss
import plot_tools as pt
import matplotlib.pyplot as plt
import json
import numpy as np
import dm_tools as dm
from read_bands import read_fort25
from copy import deepcopy
from imp import reload
pt.matplotlib_header()
ev = 27.21138386

def find_gap(bands,efermi):
  val, cond = -1,-1
  bmaxs = np.array([band.max() for band in bands])
  bmins = np.array([band.min() for band in bands])
  filled = bmaxs < efermi
  empty  = bmins > efermi
  metal = filled*empty
  if metal.any():
    return 0.0
  else:
    return bmins[empty].min() - bmaxs[filled].max()

def find_gap_k(bands_at_k,efermi):
  val, cond = -1,-1
  filled = bands_at_k < efermi
  empty  = bands_at_k > efermi
  return bands_at_k[empty].min() - bands_at_k[filled].max()

def plot_bandstructure(fort25filename):
  dat = read_fort25(fort25filename)
  efermi = dat['bands'][0]['efermi']
  print("efermi",efermi)
  fcc_kname = {
      '000':r"$\Gamma$",
      '044':r"$X$",
      '444':r"$L$",
      '264':r"$W$",
      '255':r"$U$",
      '363':r"$K$"
    }
  tmp = {
      'kpt':[],
      'klabpos':[0.0],
      'klabs':[],
      'bands':[]
    }

  pdat = {'up':deepcopy(tmp),'down':deepcopy(tmp)}

  # Classify and process.
  sp = 'up'
  for bi,b in enumerate(dat['bands']):
    if bi >= len(dat['bands'])/2: sp = 'down'
    pdat[sp]['klabs'].append(fcc_kname[b['k0']])
    pdat[sp]['klabpos'].append(b['dkp']*b['dat'].shape[0])
    pdat[sp]['kpt'] += np.repeat(b['dkp'],b['dat'].shape[0]).tolist()
    pdat[sp]['bands'] += b['dat'].tolist()

  # Make plot-friendly.
  for sp in ['up','down']:
    pdat[sp]['klabs'].append(fcc_kname[dat['bands'][-1]['k1']])
    pdat[sp]['klabpos'] = np.array(pdat[sp]['klabpos']).cumsum()
    pdat[sp]['kpt'] = np.array(pdat[sp]['kpt']).cumsum()
    pdat[sp]['bands'] = np.array(pdat[sp]['bands']).T

  fig,ax = plt.subplots(1,1)
  ax.axhline(0.0,color='k',lw=1)

  for bx,band in enumerate(pdat[sp]['bands']):
    ax.plot(pdat[sp]['kpt'],(band - efermi)*ev,label=bx)
  ax.set_xticks(pdat[sp]['klabpos'])
  ax.set_xticklabels(pdat[sp]['klabs'])
  pt.fix_lims(ax)
  ax.legend(loc='upper left',bbox_to_anchor=(1,1))
  ax.set_ylabel('Energy (eV)')
  print("Trans. gap",find_gap(pdat[sp]['bands'],efermi)*ev)
  print("Optical gap at gamma",find_gap_k(pdat[sp]['bands'].T[0],efermi)*ev)
  return fig,ax

def import_hf_data(filename):
  """ Import the data into a dictionary. """
  dat = dm.read_ekt(open(filename,'r'))
  ekt = {}
  ekt['rho'] = dict(zip(['up','dn'],[np.array(dat['o%s'%s]) for s in ['u','d']])) 
  ekt['rer'] = dict(zip(['up','dn'],[np.array(dat['o%se'%s]) for s in ['u','d']])) 
  ekt['val'] = dict(zip(['up','dn'],[np.array(dat['v%s'%s]) for s in ['u','d']])) 
  ekt['ver'] = dict(zip(['up','dn'],[np.array(dat['v%se'%s]) for s in ['u','d']])) 
  ekt['con'] = dict(zip(['up','dn'],[np.array(dat['c%s'%s]) for s in ['u','d']])) 
  ekt['cer'] = dict(zip(['up','dn'],[np.array(dat['c%se'%s]) for s in ['u','d']])) 
  return ekt

def import_data(filename):
  """ Import the data into a dictionary. """
  dat = json.load(open(filename,'r'))
  ekt = {}
  ekt['rho'] = dict(zip(['up','dn'],[np.array(dat['properties']['EKT']['obdm'][s]) for s in ['up','down']]))
  ekt['rer'] = dict(zip(['up','dn'],[np.array(dat['properties']['EKT']['obdm'][s]) for s in ['up_err','down_err']]))
  ekt['val'] = dict(zip(['up','dn'],[np.array(dat['properties']['EKT']['valence'][s]) for s in ['up','down']]))
  ekt['ver'] = dict(zip(['up','dn'],[np.array(dat['properties']['EKT']['valence'][s]) for s in ['up_err','down_err']]))
  ekt['con'] = dict(zip(['up','dn'],[np.array(dat['properties']['EKT']['conduction'][s]) for s in ['up','down']]))
  ekt['cer'] = dict(zip(['up','dn'],[np.array(dat['properties']['EKT']['conduction'][s]) for s in ['up_err','down_err']]))
  return ekt

def compute_spectra(ekt,nstats=1000):
  """ Compute the ionization and electron affinity spectra from above data. """
  norb = ekt['norb']

  boots = my.Bootstrapper_eigh(
      ekt['val']['up'][:norb,:norb],
      ekt['ver']['up'][:norb,:norb],
      ekt['rho']['up'][:norb,:norb],
      ekt['rer']['up'][:norb,:norb])
  ionboots = boots.gen_stats(nstats)
  ekt['ion'] = sci.linalg.eigh(
      ekt['val']['up'][:norb,:norb],
      ekt['rho']['up'][:norb,:norb]
    )[0]
  ekt['ier'] = ionboots['variance']**0.5

  boots = my.Bootstrapper_eig(
      ekt['con']['up'][norb:2*norb,norb:2*norb],
      ekt['cer']['up'][norb:2*norb,norb:2*norb],
      1.-ekt['rho']['up'][norb:2*norb,norb:2*norb],
      1.-ekt['rer']['up'][norb:2*norb,norb:2*norb])
  affboots = boots.gen_stats(nstats)
  ekt['aff'] = sci.linalg.eig(
      ekt['con']['up'][norb:2*norb,norb:2*norb],
      1.-ekt['rho']['up'][norb:2*norb,norb:2*norb]
    )[0]
  print(ekt['aff'])
  ekt['aer'] = affboots['variance']**0.5

  return ekt

def plot_ekt(ekt):
  fig,ax = plt.subplots(2,2,sharex=True,sharey=False)
  norb = ekt['norb']
  # Diagonal approximation.
  ax[0,0].errorbar(
      range(ekt['norb']),
      ekt['val']['up'].diagonal()[:norb]*ev,
      ekt['ver']['up'].diagonal()[:norb]*ev,
      fmt='none',ecolor=pt.pc['b'],
      capthick=1,capsize=2)
  ax[0,0].plot(
      range(ekt['norb']),
      ekt['val']['up'].diagonal()[:norb]*ev,
      'o',color='none',
      markeredgewidth=1,
      markeredgecolor=pt.pc['b'],
      label='Val.')
  ax[0,0].errorbar(
      range(ekt['norb']),
      ekt['con']['up'].diagonal()[norb:2*norb]*ev,
      ekt['cer']['up'].diagonal()[norb:2*norb]*ev,
      fmt='none',ecolor=pt.pc['r'],
      capthick=1,capsize=2)
  ax[0,0].plot(
      range(ekt['norb']),
      ekt['con']['up'].diagonal()[norb:2*norb]*ev,
      'o',color='none',
      markeredgewidth=1,
      markeredgecolor=pt.pc['r'],
      label='Cond.')
  ax[0,1].errorbar(
      range(ekt['norb']),
      ekt['con']['up'].diagonal()[norb:2*norb]*ev - \
      ekt['val']['up'].diagonal()[:norb]*ev,
      (ekt['ver']['up'].diagonal()[:norb]**2 + \
      ekt['cer']['up'].diagonal()[norb:2*norb]**2)**0.5*ev,
      fmt='none',ecolor=pt.pc['p'],
      capthick=1,capsize=2)
  ax[0,1].plot(
        range(ekt['norb']),
        ekt['con']['up'].diagonal()[norb:2*norb]*ev - \
        ekt['val']['up'].diagonal()[:norb]*ev,
        'o',color='none',
        markeredgewidth=1,
        markeredgecolor=pt.pc['p'],
        label='Cond.')

  # Full EKT.
  ax[1,0].errorbar(
      range(ekt['norb']),
      ekt['ion']*ev,
      ekt['ier']*ev,
      fmt='none',ecolor=pt.pc['b'],
      capthick=1,capsize=2)
  ax[1,0].plot(
      range(ekt['norb']),
      ekt['ion']*ev,
      'o',color='none',
      markeredgewidth=1,
      markeredgecolor=pt.pc['b'])
  ax[1,0].errorbar(
      range(ekt['norb']),
      ekt['aff'].real*ev,
      ekt['aer'].real*ev,
      fmt='none',ecolor=pt.pc['r'],
      capthick=1,capsize=2)
  ax[1,0].plot(
      range(ekt['norb']),
      ekt['aff'].real*ev,
      'o',color='none',
      markeredgewidth=1,
      markeredgecolor=pt.pc['r'])
  ax[1,1].errorbar(
      range(ekt['norb']),
      (ekt['aff'].real+ekt['ion'])*ev,
      (ekt['aer'].real**2+ekt['ier']**2)**0.5*ev,
      fmt='none',ecolor=pt.pc['p'],
      capthick=1,capsize=2)
  ax[1,1].plot(
      range(ekt['norb']),
      (ekt['aff'].real+ekt['ion'])*ev,
      'o',color='none',
      markeredgewidth=1,
      markeredgecolor=pt.pc['p'])
  for a in ax.T:
    a[-1].set_xlabel("Value index")
  for a in ax.flatten():
    pt.fix_lims(a)
  ax[0,0].set_ylabel("Matrix element (eV)")
  ax[1,0].set_ylabel("Eigenvalue (eV)")

  fig.set_size_inches(6,5)
  ax[0,0].set_title("Total values")
  ax[0,1].set_title("Sum")
  ax[0,0].legend(loc='best',frameon=True)

  return fig,ax
