#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""This modules manages the I/O to access the data"""
# This Module should be moved to snfdata when avialable

import os
import numpy as np
import pandas as pd

__all__ = ["get_source", "analysis_dataframe"]


# TMP
DATAPATH = os.path.dirname(__file__)+"/data/"
PRODNAME = "ALLEG"
#

def analysis_dataframe(sources=["hubblizer","idr","host","phrenoatmax"],
                       source_index="hubblizer"):
    """
    """
    keys_to_load = {
    "idr":       ["host.zcmb","host.zhelio",'host.zhelio.err',
                  "idr.saltprefix", "idr.subset",
                  "salt2.Color","salt2.Color.err","salt2.X1","salt2.X1.err",
                  "salt2.DayMax","salt2.DayMax.err",
                  "target.ra","target.dec","target.name","target.mwebv","target.type"],
    "hubblizer": ["hubblizer.dmfit_corr","hubblizer.dmfit_corr.err",
                  "hubblizer.dmfit_orig","hubblizer.dmfit_orig.err"],
    "phrenoatmax": ['phrenology.EWCaIIHK','phrenology.EWCaIIHK.err',
                    'phrenology.EWCaIIIR','phrenology.EWCaIIIR.err',
                    'phrenology.EWFe4800','phrenology.EWFe4800.err',
                    'phrenology.EWMgII','phrenology.EWMgII.err',
                    'phrenology.EWOI7773','phrenology.EWOI7773.err',
                    'phrenology.EWSIIW','phrenology.EWSIIW.err',
                    'phrenology.EWSIIW_L','phrenology.EWSIIW_L.err',
                    'phrenology.EWSIIW_R','phrenology.EWSIIW_R.err',
                    'phrenology.EWSiII4000','phrenology.EWSiII4000.err',
                    'phrenology.EWSiII5972','phrenology.EWSiII5972.err',
                    'phrenology.EWSiII6355','phrenology.EWSiII6355.err',
                    'phrenology.vSiII_4128','phrenology.vSiII_4128.err',
                    'phrenology.vSiII_5454','phrenology.vSiII_5454.err',
                    'phrenology.vSiII_5640','phrenology.vSiII_5640.err',
                    'phrenology.vSiII_5972','phrenology.vSiII_5972.err',
                    'phrenology.vSiII_6355','phrenology.vSiII_6355.err',
                    'phrenology.RCa','phrenology.RCa.err',
                    'phrenology.RCaIIHK','phrenology.RCaIIHK.err',
                    'phrenology.RCaIIIR','phrenology.RCaIIIR.err',
                    'phrenology.RCaS','phrenology.RCaS.err',
                    'phrenology.RCaS2','phrenology.RCaS2.err',
                    'phrenology.RFe4800','phrenology.RFe4800.err',
                    'phrenology.RMgII','phrenology.RMgII.err',
                    'phrenology.ROI7773','phrenology.ROI7773.err',
                    'phrenology.RSIIW','phrenology.RSIIW.err',
                    'phrenology.RSIIW_L','phrenology.RSIIW_L.err',
                    'phrenology.RSIIW_R','phrenology.RSIIW_R.err',
                    'phrenology.RSi','phrenology.RSi.err',
                    'phrenology.RSiII4000','phrenology.RSiII4000.err',
                    'phrenology.RSiII5972','phrenology.RSiII5972.err',
                    'phrenology.RSiII6355','phrenology.RSiII6355.err',
                    'phrenology.RSiS','phrenology.RSiS.err',
                    'phrenology.RSiSS','phrenology.RSiSS.err',
                    'phrenology.Rsjb','phrenology.Rsjb.err'],
        "host": ["localmass","localmass.err","sfr","sfr.err",
                 "dispersion","dispersion.err",
                 "HA","HA.err","NII","NII.err","OII","OII.err",'localradius_kpc',
                 'uly_age', 'uly_chi2','uly_feh',
                 'velocity', 'velocity.err', u'zhelio', u'zhelio.err',
                 'zhelio.quality', 'zhelio.quality.info']
        }
        
    ## Load the sources and the analysis_frame
    sourceframe = {}
    for sourcename in sources:
        sourceframe[sourcename] = pd.DataFrame(get_source(sourcename)).T

    analysis_frame = pd.DataFrame(index=sourceframe[source_index].index)

    for sourcename in sources:
        analysis_frame[keys_to_load[sourcename]] = sourceframe[sourcename][keys_to_load[sourcename]]
    
    return analysis_frame
        
    
# ----------------------- #
# - Individual GET      - #
# ----------------------- #
def get_source(sourcename,**kwargs):
    """ source name must be a known snfanalysis source:
       idr/hubblizer/host/phreno/phrenoatmax
    """
    if sourcename.lower() in ["idr","meta"]:
        return get_snidr()
    if sourcename.lower() in ["hubblizer","hr"]:
        return get_snhubblizer()
    if sourcename.lower() in ["host","localhost","idrhost"]:
        return get_hostidr()
    if sourcename.lower() in ["phreno","phrenology"]:
        return get_phreno(atmax=False)
    if sourcename.lower() in ["phrenoatmax"]:
        return get_phreno(atmax=True,raw=False)
    raise ValueError("sourcename %s not recongnized (use: idr/hubblizer/host/phreno/phrenoatmax)"%sourcename)

def get_snidr():
    """ IDR of the sn data, including SALT """
    return load_pkl(DATAPATH+PRODNAME+"/META.pkl")

def get_snhubblizer():
    """ IDR of the sn data, including SALT and Hubble Residuals """
    return load_pkl(DATAPATH+PRODNAME+"/SNF-0203-ALLEG2a_SNeIa_hubble.pkl")

def get_hostidr():
    """ IDR containing the host information """
    return load_pkl(DATAPATH+PRODNAME+"/localhost_idr.pkl")

def get_phreno(atmax=False, raw=True):
    if atmax:
        atmaxphreno = load_pkl(DATAPATH+PRODNAME+"/phreno_ALLEG2_20160414__plusorminus2.5daysatmax_daymaxerr_lt_1.pkl")
        if raw:
            return atmaxphreno
        atmaxphreno_flat = {}
        for name, d in atmaxphreno.items():
            if len(d.keys())<2:
                atmaxphreno_flat[name] = d.values()[0]
            else:
                minphase = np.argmin([np.abs(dd["salt2.phase"]) for dd in d.values()])
                atmaxphreno_flat[name] = d.values()[minphase]
        return atmaxphreno_flat
                 
        
    return load_pkl(DATAPATH+PRODNAME+"/phreno_ALLEG2_20160414.pkl")





# ----------------------- #
# - ToolBox             - #
# ----------------------- #
def load_pkl(filename):
    """
    """
    import cPickle as pkl
    try:
        pkl_file = open(filename,'rb')
    except:
        raise IOError("The given file does not exist %s"%filename)
    
    return pkl.load(pkl_file)


def dump_pkl(data,filename,**kwargs):
    """
    """
    from cPickle import dump
    if len(filename.split("."))>1 and filename.split(".")[-1]=="pkl":
        outfile =  open(filename,"wb")
    else:
        outfile =  open(filename+".pkl","wb")
    
    dump(data, outfile,**kwargs)
    outfile.close()
