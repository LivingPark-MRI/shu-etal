import SimpleITK as sitk
import numpy, math
from radiomics.glcm import RadiomicsGLCM
from radiomics.glrlm import RadiomicsGLRLM
from radiomics.glszm import RadiomicsGLSZM

'''
Required to save 'angles' as a variable in the GLCM & GLRLM classes
'''

class RadiomicsHelper:
    '''
    Parameters
    ----------
    image : sitk.Image
        MRI image of patient
    mask : sitk.Image
        White matter mask of patient
    distance : int
        Offset to calculate matrix, single number, must be a list
    binCount: int
        Binning level to use for pre-processing
    '''
    def __init__(self, image, mask, distance=[1], binCount=32):
        if not isinstance(image, sitk.Image):
            raise ValueError("Image must be a SimpleITK (sitk) Image.")
        if not isinstance(mask, sitk.Image):
            raise ValueError("Mask must be a SimpleITK (sitk) Image.")
        if not isinstance(distance, list):
            raise ValueError("Distance must be a list.")
        if len(distance) != 1:
            raise ValueError("Only a single distance can be defined.")
        if not isinstance(binCount, int):
            raise ValueError("BinCount must be an integer.")

        self.image = image
        self.mask = mask
        self.distance = distance
        self.binCount = binCount

        # Initialize GLCM, GLRLM, GLSZM matrixes
        try:
            self.glcm = RadiomicsGLCM(image, mask, binCount=binCount, distances=distance)
            self.glrlm = RadiomicsGLRLM(image, mask, binCount=binCount, distances=distance)
            self.glszm = RadiomicsGLSZM(image, mask, binCount=binCount, distances=distance)
            
            self.glcm._initCalculation()
            self.glrlm._initCalculation()
            self.glszm._initCalculation()

            if not hasattr(self.glcm, 'angles'):
                raise ValueError("Error, need to save angles as a variable in GLCM.py")
            if not hasattr(self.glrlm, 'angles'):
                raise ValueError("Error, need to save angles as a variable in GLRLM.py")

            if not self.verifyValues():
                raise Exception("Error found in test cases, please verify.")
            
        except Exception as e:
            print(f"Error initializing glcm: {e}")
    
    '''
    =============
    GLCM features
    =============
    '''

    def getClusterShadeFeatureValue(self, angle):
        i = self.glcm.coefficients['i']
        j = self.glcm.coefficients['j']
        ux = self.glcm.coefficients['ux'][:,:,:,angle]
        uy = self.glcm.coefficients['uy'][:,:,:,angle]
        cs = numpy.sum((self.glcm.P_glcm[:, :, :, angle] * (((i + j)[None, :, :] - ux - uy) ** 3)), (1, 2))
        return numpy.nanmean(cs)

    def getClusterProminenceFeatureValue(self, angle):
        i = self.glcm.coefficients['i']
        j = self.glcm.coefficients['j']
        ux = self.glcm.coefficients['ux'][:,:,:,angle]
        uy = self.glcm.coefficients['uy'][:,:,:,angle]
        cp = numpy.sum((self.glcm.P_glcm[:, :, :, angle] * (((i + j)[None, :, :] - ux - uy) ** 4)), (1, 2))
        return numpy.nanmean(cp)

    def getContrastFeatureValue(self, angle):
        i = self.glcm.coefficients['i']
        j = self.glcm.coefficients['j']
        cont = numpy.sum((self.glcm.P_glcm[:, :, :, angle] * ((numpy.abs(i - j))[None, :, :] ** 2)), (1, 2))
        return numpy.nanmean(cont)            

    def getAutocorrelationFeatureValue(self, angle):
        i = self.glcm.coefficients['i']
        j = self.glcm.coefficients['j']
        ac = numpy.sum(self.glcm.P_glcm[:, :, :, angle] * (i * j)[None, :, :], (1, 2))
        return numpy.nanmean(ac)

    def getCorrelationFeatureValue(self, angle):
        eps = self.glcm.coefficients['eps']
        i = self.glcm.coefficients['i']
        j = self.glcm.coefficients['j']
        ux = self.glcm.coefficients['ux']
        uy = self.glcm.coefficients['uy']
    
        # shape = (Nv, 1, 1, angles)
        sigx = numpy.sum(self.glcm.P_glcm[:, :, :, angle] * ((i[None, :, :, None] - ux) ** 2)[:,:,:,angle], (1, 2), keepdims=True) ** 0.5
        # shape = (Nv, 1, 1, angles)
        sigy = numpy.sum(self.glcm.P_glcm[:, :, :, angle] * ((j[None, :, :, None] - uy) ** 2)[:,:,:,angle], (1, 2), keepdims=True) ** 0.5
    
        corm = numpy.sum(self.glcm.P_glcm[:, :, :, angle] * (i[None, :, :, None] - ux)[:,:,:,angle] * (j[None, :, :, None] - uy)[:,:,:,angle], (1, 2), keepdims=True)
        corr = corm / (sigx * sigy + eps)
        corr[sigx * sigy == 0] = 1  # Set elements that would be divided by 0 to 1.
        return numpy.nanmean(corr)
    
    '''
    =============
    GLRLM features
    =============
    '''

    # def getHighGrayLevelRunEmphasisFeatureValue(self, angle):
    #     pg = self.glrlm.coefficients['pg'][:, :, angle]
    #     ivector = self.glrlm.coefficients['ivector']
    #     Nr = self.glrlm.coefficients['Nr'][:, angle]

    #     hglre = numpy.sum((pg * (ivector[None, :, None] ** 2)), 1) / Nr
    #     return numpy.nanmean(hglre)

    # def getLongRunEmphasisFeatureValue(self, angle):
    #     pr = self.glrlm.coefficients['pr'][:, :, angle]
    #     jvector = self.glrlm.coefficients['jvector']
    #     Nr = self.glrlm.coefficients['Nr'][:, angle]

    #     lre = numpy.sum((pr * (jvector[None, :, None] ** 2)), 1) / Nr
    #     return numpy.nanmean(lre)

    # def getLongRunHighGrayLevelEmphasisFeatureValue(self, angle):
    #     ivector = self.glrlm.coefficients['ivector']
    #     jvector = self.glrlm.coefficients['jvector']
    #     Nr = self.glrlm.coefficients['Nr'][:, angle]

    #     lrhgle = numpy.sum((self.glrlm.P_glrlm * ((jvector[None, None, :, None] ** 2) * (ivector[None, :, None, None] ** 2))),
    #                     (1, 2)) / Nr
    #     return numpy.nanmean(lrhgle)   

    # def getShortRunEmphasisFeatureValue(self, angle):
    #     pr = self.glrlm.coefficients['pr'][:, :, angle]
    #     jvector = self.glrlm.coefficients['jvector']
    #     Nr = self.glrlm.coefficients['Nr'][:, angle]

    #     sre = numpy.sum((pr / (jvector[None, :, None] ** 2)), 1) / Nr
    #     return numpy.nanmean(sre)

    # def getShortRunHighGrayLevelEmphasisFeatureValue(self, angle):
    #     ivector = self.glrlm.coefficients['ivector']
    #     jvector = self.glrlm.coefficients['jvector']
    #     Nr = self.glrlm.coefficients['Nr'][:, angle]

    #     srhgle = numpy.sum((self.glrlm.P_glrlm * (ivector[None, :, None, None] ** 2) / (jvector[None, None, :, None] ** 2)),
    #                     (1, 2)) / Nr
    #     return numpy.nanmean(srhgle)

    def getLongRunLowGrayLevelEmphasisFeatureValue(self, angle):
        ivector = self.glrlm.coefficients['ivector']
        jvector = self.glrlm.coefficients['jvector']
        Nr = self.glrlm.coefficients['Nr'][:, angle]

        lrlgle = numpy.sum((self.glrlm.P_glrlm * (jvector[None, None, :, None] ** 2) / (ivector[None, :, None, None] ** 2)),
                        (1, 2)) / Nr
        return numpy.nanmean(lrlgle)

    def getGrayLevelNonUniformityFeatureValue(self, angle):
        pg = self.glrlm.coefficients['pg'][:, :, angle]
        Nr = self.glrlm.coefficients['Nr'][:, angle]

        gln = numpy.sum((pg ** 2), 1) / Nr
        return numpy.nanmean(gln)  

    def getRunLengthNonUniformityFeatureValue(self, angle):
        pr = self.glrlm.coefficients['pr'][:, :, angle]
        Nr = self.glrlm.coefficients['Nr'][:, angle]
    
        rln = numpy.sum((pr ** 2), 1) / Nr
        return numpy.nanmean(rln)

    def getShortRunLowGrayLevelEmphasisFeatureValue(self, angle):
        ivector = self.glrlm.coefficients['ivector']
        jvector = self.glrlm.coefficients['jvector']
        Nr = self.glrlm.coefficients['Nr'][:, angle]

        srlgle = numpy.sum((self.glrlm.P_glrlm[:, :, :, angle] / ((ivector[None, :, None] ** 2) * (jvector[None, None, :] ** 2))),
                           (1, 2)) / Nr
        return numpy.nanmean(srlgle)
    
    '''
    =============
    General utils
    =============
    '''

    def getAngleIndex(self, isGLCM=True):
        '''
        Returns dictionary of 2D angle indexes
        '''
        angleDict= {}
        angles = self.glcm.angles if isGLCM else self.glrlm.angles

        def convertToDegree(x, y):
            if y != 0:
                return math.degrees(math.atan(-x/y))%180
            else:
                return float(90%180)

        for index, angle in enumerate(angles):
            x, y, z = angle
            
            if z==0:
                degree = convertToDegree(x,y)
                if degree in [float(0),float(45), float(90), float(135)]:
                    angleDict[int(degree)] = index

        return angleDict

    def getAllFeatures(self):
        features = {}
        glcmAngles = self.getAngleIndex()
        glrlmAngles = self.getAngleIndex(isGLCM=False)
        distance = self.distance[0]

        for angle in glcmAngles:
            features[f"Inertia_angle{angle}_offset{distance}"] = self.getContrastFeatureValue(glcmAngles[angle])
            features[f"HaralickCorrelation_angle{angle}_offset{distance}"] = self.getCorrelationFeatureValue(glcmAngles[angle])
            features[f"ClusterProminenceFeatureValue_angle{angle}_offset{distance}"] = self.getClusterProminenceFeatureValue(glcmAngles[angle])
            features[f"ClusterShadeFeatureValue_angle{angle}_offset{distance}"] = self.getClusterShadeFeatureValue(glcmAngles[angle])
            features[f"DifferenceAverageFeatureValue_AllDirection_offset{distance}"] = self.glcm.getDifferenceAverageFeatureValue()[0]
            features[f"GLCMEntropy_AllDirection_offset{distance}"] = self.glcm.getJointEntropyFeatureValue()[0]

        for angle in glrlmAngles:
            features[f"LongRunLowGrayLevelEmphasis_angle{angle}_offset{distance}"] = self.getLongRunLowGrayLevelEmphasisFeatureValue(glrlmAngles[angle])
            features[f"GrayLevelNonUniformity_angle{angle}_offset{distance}"] = self.getGrayLevelNonUniformityFeatureValue(glrlmAngles[angle])
            features[f"RunLengthNonUniformity_angle{angle}_offset{distance}"] = self.getRunLengthNonUniformityFeatureValue(glrlmAngles[angle])
            features[f"ShortRunLowGrayLevelEmphasis_angle{angle}_offset{distance}"] = self.getShortRunLowGrayLevelEmphasisFeatureValue(glrlmAngles[angle])

        features[f"SizeZoneNonUniformity_offset{distance}"] = self.glszm.getSizeZoneNonUniformityFeatureValue()[0]
        features[f"LargeAreaEmphasisFeatureValue_offset{distance}"] = self.glszm.getLargeAreaEmphasisFeatureValue()[0]
        features[f"SmallAreaEmphasisFeatureValue_offset{distance}"] = self.glszm.getSmallAreaEmphasisFeatureValue()[0]
        features[f"ZonePercentageFeatureValue_offset{distance}"] = self.glszm.getZonePercentageFeatureValue()[0]

        return features

    '''
    =============
    UNIT TESTING
    =============
    '''

    def verifyValues(self):
        
        def compute_average_over_angles(angles, fun):
            s = []
            for i in range(angles):
                s.append(fun(i))
            return sum(s) / len(s)

        glcmAngles = self.glcm.P_glcm.shape[-1] 
        glrmAngles = self.glrlm.P_glrlm.shape[-1]

        # Contrast | GLCM
        avg_sum = compute_average_over_angles(glcmAngles, self.getContrastFeatureValue)
        assert (abs(avg_sum - self.glcm.getContrastFeatureValue()[0]) < 0.0005), f"Contrast not equal! Average over angles is {avg_sum} while computed value is {self.glcm.getContrastFeatureValue()[0]}"

        # Correlation (inertia) | GLCM
        avg_sum = compute_average_over_angles(glcmAngles, self.getAutocorrelationFeatureValue)
        assert (abs(avg_sum - self.glcm.getAutocorrelationFeatureValue()[0]) < 0.0005), f"Autocorrelation not equal! Average over angles is {avg_sum} while computed value is {self.glcm.getAutocorrelationFeatureValue()[0]}"
        
        # getClusterShadeFeatureValue | GLCM
        avg_sum = compute_average_over_angles(glcmAngles, self.getClusterShadeFeatureValue)
        assert (abs(avg_sum - self.glcm.getClusterShadeFeatureValue()[0]) < 0.0005), f"ClusterShade not equal! Average over angles is {avg_sum} while computed value is {self.glcm.getClusterShadeFeatureValue()[0]}"
        
        # getClusterProminenceFeatureValue | GLCM
        avg_sum = compute_average_over_angles(glcmAngles, self.getClusterProminenceFeatureValue)
        assert (abs(avg_sum - self.glcm.getClusterProminenceFeatureValue()[0]) < 0.0005), f"ClusterProminence not equal! Average over angles is {avg_sum} while computed value is {self.glcm.getClusterProminenceFeatureValue()[0]}"

        # getLongRunLowGrayLevelEmphasisFeatureValue | GLRLM
        avg_sum = compute_average_over_angles(glrmAngles, self.getLongRunLowGrayLevelEmphasisFeatureValue)
        assert (abs(avg_sum - self.glrlm.getLongRunLowGrayLevelEmphasisFeatureValue()[0]) < 0.0005), f"LongRunLowGrayLevelEmphasis not equal! Average over angles is {avg_sum} while computed value is {self.glrlm.getLongRunLowGrayLevelEmphasisFeatureValue()[0]}"

        # getGrayLevelNonUniformityFeatureValue | GLRLM
        avg_sum = compute_average_over_angles(glrmAngles, self.getGrayLevelNonUniformityFeatureValue)
        assert (abs(avg_sum - self.glrlm.getGrayLevelNonUniformityFeatureValue()[0]) < 0.0005), f"GrayLevelNonUniformity not equal! Average over angles is {avg_sum} while computed value is {self.glrlm.getGrayLevelNonUniformityFeatureValue()[0]}"

        # getRunLengthNonUniformityFeatureValue | GLRLM
        avg_sum = compute_average_over_angles(glrmAngles, self.getRunLengthNonUniformityFeatureValue)
        assert (abs(avg_sum - self.glrlm.getRunLengthNonUniformityFeatureValue()[0]) < 0.0005), f"RunLength not equal! Average over angles is {avg_sum} while computed value is {self.glrlm.getRunLengthNonUniformityFeatureValue()[0]}"

        # getRunLengthNonUniformityFeatureValue
        avg_sum = compute_average_over_angles(glrmAngles, self.getShortRunLowGrayLevelEmphasisFeatureValue)
        assert (abs(avg_sum - self.glrlm.getShortRunLowGrayLevelEmphasisFeatureValue()[0]) < 0.0005), f"Short Run Low Gray Level Emphasis not equal! Average over angles is {avg_sum} while computed value is {self.glrlm.getShortRunLowGrayLevelEmphasisFeatureValue()[0]}"

        # getCorrelationFeatureValue
        avg_sum = compute_average_over_angles(glcmAngles, self.getCorrelationFeatureValue)
        assert (abs(avg_sum - self.glcm.getCorrelationFeatureValue()[0]) < 0.0005), f"Correlation not equal! Average over angles is {avg_sum} while computed value is {self.glcm.getCorrelationFeatureValue()[0]}"

        return True