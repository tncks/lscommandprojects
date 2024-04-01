#!/usr/bin/env python3
# -*- coding: utf-8 -*
# sample_python aims to allow seamless integration with lua.
# see examples below

import os
import sys
import pdb  # use pdb.set_trace() for debugging
import code # or use code.interact(local=dict(globals(), **locals()))  for debugging.
import xml.etree.ElementTree as ET
import numpy as np
import math
from PIL import Image 
class Color:
    def __init__(self, R, G, B):
        self.color=np.array([R,G,B]).astype(np.float64)

    # Gamma corrects this color.
    # @param gamma the gamma value to use (2.2 is generally used).
    def gammaCorrect(self, gamma):
        inverseGamma = 1.0 / gamma;
        self.color=np.power(self.color, inverseGamma)

    def toUINT8(self):
        return (np.clip(self.color, 0,1)*255).astype(np.uint8)



def main():


    tree = ET.parse(sys.argv[1])
    root = tree.getroot()

    # set default values
    #print(redeem['_zero'])
    viewDir=np.array([0,0,-1]).astype(np.float64)
    viewUp=np.array([0,1,0]).astype(np.float64)
    viewProjNormal=-1*viewDir  # you can safely assume this. (no examples will use shifted perspective camera)
    viewWidth=1.0
    viewHeight=1.0
    projDistance=1.0
    intensity=np.array([1,1,1]).astype(np.float64)  # how bright the light is.
    #print(np.cross(viewDir, viewUp))

    imgSize=np.array(root.findtext('image').split()).astype(np.int32)
    #print('Output expectation sized: ',imgSize[0])
    #print('with by height: ',imgSize[1])

    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    for c in root.findall('camera'):
        viewPoint=np.array(c.findtext('viewPoint').split()).astype(np.float64)
        viewDir_=np.array(c.findtext('viewDir').split()).astype(np.float64)
        viewProjNormal_=np.array(c.findtext('projNormal').split()).astype(np.float64)
        viewUp_=np.array(c.findtext('viewUp').split()).astype(np.float64)
        projDistance_=np.array(c.findtext('projDistance').split()).astype(np.float64)
        viewWidth_=np.array(c.findtext('viewWidth').split()).astype(np.float64)
        viewHeight_=np.array(c.findtext('viewHeight').split()).astype(np.float64)
        print('viewpoint', viewPoint)

    # Create an empty image
    channels=3
    img = np.zeros((imgSize[1], imgSize[0], channels), dtype=np.uint8)
    img[:,:]=0
    
    # replace the code block below!
    w = viewDir_
    u = np.cross(w, viewUp_)
    v = np.cross(w, u)

    for x in np.arange(imgSize[0]):
        for y in np.arange(imgSize[1]):
            ray = 

    rawimg = Image.fromarray(img, 'RGB')
    #rawimg.save('out.png')
    rawimg.save(sys.argv[1]+'.png')
    
if __name__=="__main__":
    main()
