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

def foo(cen, r, e, d):
    
    a = np.linalg.norm(d) ** 2  
    b = 2 * np.dot(d, e - cen)
    c = np.linalg.norm(e - cen) ** 2 - r**2
    discrim = b**2 - 4 * a * c

    if discrim > 0:
        t1 = (-b - np.sqrt(discrim)) / (2 * a)
        t2 = (-b + np.sqrt(discrim)) / (2 * a)

        if t1 > 0 and t2 > 0:  
            return min(t1, t2)  

    return None

def goo(objects, e, d, cen, r):
    
    t_distances_set = [
        foo(cen, r, e, d)
        #foo(obj["center"], obj["radius"], e, d)
        #for obj in objects
    ]
    nearest_object = None
    min_t = np.inf

    for index, distance in enumerate(t_distances_set):
        if distance and distance < min_t:
            min_t = distance
            nearest_object = True # temporal flag substituted for testing
            #nearest_object = objects[index]

    return nearest_object, min_t


def main():


    tree = ET.parse(sys.argv[1])
    root = tree.getroot()

    # set default values
    redeem = {'_zero': '_0'} # this is a dictionary data for naming space temporarily.
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
        projNormal_=np.array(c.findtext('projNormal').split()).astype(np.float64)
        viewUp_=np.array(c.findtext('viewUp').split()).astype(np.float64)
        projDistance_=np.array(c.findtext('projDistance').split())[0].astype(np.float64)
        viewWidth_=np.array(c.findtext('viewWidth').split())[0].astype(np.float64)
        viewHeight_=np.array(c.findtext('viewHeight').split())[0].astype(np.float64)
        print('viewpoint', viewPoint)
        print('viewdir', viewDir_)
        print('projNormal', projNormal_)
        print('viewUp', viewUp_)
        print('projDistance', projDistance_)
        print('viewWidth', viewWidth_)
        print('viewHeight', viewHeight_)
        redeem['e'] = viewPoint
        redeem['e_x'] = viewPoint[0]
        redeem['e_y'] = viewPoint[1]
        redeem['e_z'] = viewPoint[2]
        redeem['d'] = viewDir_
        print('projNormal', projNormal_)
        print('viewUp', viewUp_)
        print('projDistance', projDistance_)
        print('viewWidth', viewWidth_)
        print('viewHeight', viewHeight_)

        #print('redeem all inter:', redeem)

    for c in root.findall('shader'):
        diffuseColor_c=np.array(c.findtext('diffuseColor').split()).astype(np.float64)
        specularColor_=np.array(c.findtext('specularColor').split()).astype(np.float64)
        exponent_=np.array(c.findtext('exponent').split())[0].astype(np.uint8)
        #print('name', c.get('name'))
        #print('type', c.get('type'))
        print('diffuseColor', diffuseColor_c)
        print('specularColor', specularColor_)
        print('exponent', exponent_)
    for c in root.findall('surface'):
        center_=np.array(c.findtext('center').split()).astype(np.float64)
        radius_=np.array(c.findtext('radius').split())[0].astype(np.float64)
        print('type', c.get('type'))
        print('center', center_)
        print('radius', radius_)
        redeem['target_type'] = c.get('type')
        redeem['c'] = center_
        redeem['r'] = radius_
    for c in root.findall('light'):
        position_=np.array(c.findtext('position').split()).astype(np.float64)
        print('position', position_)
        intensity_=np.array(c.findtext('intensity').split()).astype(np.float64)
        print('intensity', intensity_)
        redeem['position'] = position_
        #redeem['r'] = radius_

        print('redeem all inter:', redeem)
    #code.interact(local=dict(globals(), **locals()))  
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    

    # Create an empty image
    channels=3
    img = np.zeros((imgSize[1], imgSize[0], channels), dtype=np.uint8)
    img[:,:]=0
    
    # replace the code block below!
    blue = Color(0,0,1)
    '''
    for i in np.arange(imgSize[1]): 
        white=Color(1,1,1)
        red=Color(1,0,0)
        blue=Color(0,0,1)
        img[10][i]=white.toUINT8()
        img[i][i]=red.toUINT8()
        img[i][0]=blue.toUINT8()

    for x in np.arange(imgSize[0]): 
        img[5][x]=[255,255,255]
    '''
    for i in np.arange(imgSize[1]):
        for j in np.arange(imgSize[0]):
            pixel = np.array([j,i,0])
            d = pixel - redeem['e']
            d /= np.linalg.norm(d)

            nearest_object_dummy_kbs, min_t = goo([],redeem['e'], d, redeem['c'], redeem['r'])

            if nearest_object_dummy_kbs is None:
                continue


            intersection = redeem['e'] + min_t * d

            #added line for temporal safety check
            if math.isinf(0-intersection[0]):
                #print('warning: intersection now is minus inf so continue would be!')
                continue
            #added line for temporal safety check
            if math.isinf(0-intersection[1]):
                #print('warning: intersection now is minus inf so continue would be!')
                continue
            #added line for temporal safety check
            if math.isinf(0-intersection[2]):
                #print('warning: intersection now is minus inf so continue would be!')
                continue

            print('[[IMPORTANT]] intersection: ', intersection)
            print('min_t: ', min_t)
            print('d: ', d)
            #print('LASTEOF]].. redeemed in D value: ', redeem['d'])
            #intersection = redeem['e'] + min_t * redeem['d']

            normal_to_surface = intersection - redeem['c']
            #normal_to_surface = intersection - nearest_object_dummy_kbs['c']
            normal_to_surface /= np.linalg.norm(normal_to_surface)

            light_direction = redeem['position'] - intersection
            light_direction /= np.linalg.norm(light_direction)

            _, min_t = goo([], intersection, light_direction, redeem['c'], redeem['r'])
            
            intersect_to_light_dist = np.linalg.norm(redeem['position'] - intersection)
            #is_shadowed = min_t < intersect_to_light_dist

            #if is_shadowed:
            #    continue

            print('intersection found counting')
            img[i][j]=blue.toUINT8()



    rawimg = Image.fromarray(img, 'RGB')
    #rawimg.save('out.png')
    rawimg.save(sys.argv[1]+'.png')
    
if __name__=="__main__":
    main()
