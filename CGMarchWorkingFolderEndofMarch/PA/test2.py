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
        self.color = np.array([R, G, B]).astype(np.float64)

    def toUINT8(self):
        return (np.clip(self.color, 0, 1) * 255).astype(np.uint8)

def ray_sphere_intersection(center, radius, origin, direction):
    oc = origin - center
    a = np.linalg.norm(direction) ** 2
    b = 2 * np.dot(oc, direction)
    c = np.linalg.norm(oc) ** 2 - radius ** 2
    discriminant = b ** 2 - 4 * a * c

    if discriminant < 0:
        return None, None

    t1 = (-b - math.sqrt(discriminant)) / (2 * a)
    t2 = (-b + math.sqrt(discriminant)) / (2 * a)

    if t1 > 0 and t2 > 0:
        intersection1 = origin + t1 * direction
        intersection2 = origin + t2 * direction
        return min(t1, t2), intersection1 if t1 < t2 else intersection2

    return None, None

def main():
    tree = ET.parse(sys.argv[1])
    root = tree.getroot()

    imgSize = np.array(root.findtext('image').split()).astype(np.int32)
    img = np.zeros((imgSize[1], imgSize[0], 3), dtype=np.uint8)

    for c in root.findall('camera'):
        viewPoint = np.array(c.findtext('viewPoint').split()).astype(np.float64)
        viewDir = np.array(c.findtext('viewDir').split()).astype(np.float64)

    for c in root.findall('light'):
        lightPos = np.array(c.findtext('position').split()).astype(np.float64)
        lightIntensity = np.array(c.findtext('intensity').split()).astype(np.float64)

    for c in root.findall('surface'):
        surfaceType = c.get('type')
        center = np.array(c.findtext('center').split()).astype(np.float64)
        radius = float(c.findtext('radius'))

        shaderRef = c.find('shader').get('ref')
        for shader in root.findall('shader'):
            if shader.get('name') == shaderRef:
                diffuseColor = np.array(shader.findtext('diffuseColor').split()).astype(np.float64)
                specularColor = np.array(shader.findtext('specularColor').split()).astype(np.float64)
                exponent = float(shader.findtext('exponent'))

        for i in range(imgSize[1]):
            for j in range(imgSize[0]):
                pixel = np.array([j, i, 0])
                rayDirection = pixel - viewPoint
                rayDirection /= np.linalg.norm(rayDirection)
                t, intersection = ray_sphere_intersection(center, radius, viewPoint, rayDirection)

                if t is None:
                    continue
                normal = intersection - center
                normal /= np.linalg.norm(normal)
                lightDir = lightPos - intersection
                lightDir /= np.linalg.norm(lightDir)
                _, min_distance = ray_sphere_intersection(center, radius, intersection, lightDir)
                intersect_to_light_dist = np.linalg.norm(lightPos - intersection)
                #is_shadowed = min_distance < intersect_to_light_dist
                is_shadowed = False
                if not is_shadowed:
                    diffuseTerm = max(0, np.dot(normal, lightDir))
                    viewDir = -rayDirection
                    reflectDir = 2 * np.dot(lightDir, normal) * normal - lightDir
                    specularTerm = max(0, np.dot(viewDir, reflectDir)) ** exponent

                    finalColor = (diffuseColor * lightIntensity * diffuseTerm + specularColor * lightIntensity * specularTerm)
                    img[i][j] = (finalColor * 255).astype(np.uint8)



    rawimg = Image.fromarray(img, 'RGB')
    rawimg.save(sys.argv[1] + '.png')

if __name__ == "__main__":
    main()

#end
