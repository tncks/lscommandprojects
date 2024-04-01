import os
import sys
import xml.etree.ElementTree as ET
import numpy as np
from PIL import Image
import math

class Color:
    def __init__(self, R, G, B):
        self.color = np.array([R, G, B]).astype(np.float64)

    # Gamma corrects this color.
    # @param gamma the gamma value to use (2.2 is generally used).
    def gammaCorrect(self, gamma):
        inverseGamma = 1.0 / gamma
        self.color = np.power(self.color, inverseGamma)

    def toUINT8(self):
        return (np.clip(self.color, 0, 1) * 255).astype(np.uint8)

def ray_sphere_intersection(center, radius, origin, direction):
    oc = origin - center
    a = np.dot(direction, direction)
    b = 2 * np.dot(oc, direction)
    c = np.dot(oc, oc) - radius * radius
    discriminant = b * b - 4 * a * c

    if discriminant < 0:
        return None

    t = (-b - math.sqrt(discriminant)) / (2 * a)

    if t > 0:
        return t

    return None

def main():
    tree = ET.parse(sys.argv[1])
    root = tree.getroot()

    imgSize = np.array(root.findtext('image').split()).astype(np.int32)
    img = np.zeros((imgSize[1], imgSize[0], 3), dtype=np.uint8)

    for c in root.findall('camera'):
        viewPoint = np.array(c.findtext('viewPoint').split()).astype(np.float64)
        viewDir = np.array(c.findtext('viewDir').split()).astype(np.float64)
        viewUp = np.array(c.findtext('viewUp').split()).astype(np.float64)
        viewWidth = float(c.findtext('viewWidth'))
        viewHeight = float(c.findtext('viewHeight'))

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

                t = ray_sphere_intersection(center, radius, viewPoint, rayDirection)
                if t is not None:
                    # Phong shading
                    intersection = viewPoint + t * rayDirection
                    normal = intersection - center
                    normal /= np.linalg.norm(normal)
                    lightDir = lightPos - intersection
                    lightDir /= np.linalg.norm(lightDir)
                    diffuseTerm = max(0, np.dot(normal, lightDir))
                    viewDir = -rayDirection
                    reflectDir = 2 * np.dot(lightDir, normal) * normal - lightDir
                    specularTerm = max(0, np.dot(viewDir, reflectDir)) ** exponent

                    finalColor = (diffuseColor * lightIntensity * diffuseTerm +
                                  specularColor * lightIntensity * specularTerm)
                    img[i][j] = (finalColor * 255).astype(np.uint8)

    rawimg = Image.fromarray(img, 'RGB')
    rawimg.save(sys.argv[1] + '.png')

if __name__ == "__main__":
    main()

