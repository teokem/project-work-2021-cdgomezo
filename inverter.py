from numpy import *
from numpy.linalg import *
from copy import *

class Inverter:
    def __init__(self,
                 trueFluxes,
                 trueInit,
                 model,
                 presFluxKey = [],
                 realObs = False,
                 fromTrue = {},
                 trueObs = None,
                 realError = False,
                 errorObs = {},
                 errorInit = {},
                 errorFluxes = {}):
        """
        Class that contains the functions necessary to convert from flux to vector and viceversa,
        observations to vector and viceversa, \Delta^14CO_2 to F^14CO2 and viceversa, and the inversion:
        Kernel matrix, Gain matrix, Cost function, and the posterior computation.
        """
        
        self.trueFluxes = trueFluxes
        self.trueInit = trueInit
        self.model = model
        self.presFluxKey = presFluxKey
        self.realObs = realObs
        self.fromTrue = fromTrue
        self.trueObs = trueObs
        self.realError = realError
        self.errorObs = errorObs
        self.errorInit = errorInit
        self.errorFluxes = errorFluxes
       
        
    def genPrior(self): 
        """
        Generates prior fluxes from true fluxes and returns the control vector.
        """
        self.priorFluxes = deepcopy(self.trueFluxes)
        self.trueVector = self.f2v(self.trueInit, self.trueFluxes)
        
        if self.realObs:
            self.stateVector = self.f2v(self.trueInit, self.priorFluxes)
        else:
            for key_i, value_i in self.priorFluxes.items():
                if key_i == 'time' or key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.priorFluxes[key_i].items():
                        for key_k, value_k in self.priorFluxes[key_i][key_j].items():
                            for key_l, value_l in self.presFluxKey.items():
                                if key_k in self.presFluxKey[key_l]:
                                    self.priorFluxes[key_i][key_j][key_k] = value_k * 0
                                else:
                                    self.priorFluxes[key_i][key_j][key_k] = value_k * self.fromTrue[key_i][key_j][key_k]

            self.stateVector = self.f2v(self.trueInit, self.priorFluxes) 

        return self.priorFluxes, self.stateVector, self.trueVector
    
    
    def genPrescribed(self):
        """
        Generates prescribed fluxes. The initial condition is prescribed by default.
        """
        prescribedInit = deepcopy(self.trueInit)

        for key_i, value_i in prescribedInit.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in prescribedInit[key_i].items():
                    prescribedInit[key_i][key_j] = value_j * 0

        self.prescribedFlux = deepcopy(self.trueFluxes)

        for key_i, value_i in self.prescribedFlux.items():
            if key_i == 'time' or key_i == 'units':
                pass
            else:
                for key_j, value_j in self.prescribedFlux[key_i].items():
                    for key_k, value_k in self.prescribedFlux[key_i][key_j].items():
                        for key_l, value_l in self.presFluxKey.items():
                            if key_k in self.presFluxKey[key_l]:
                                pass
                            else:
                                self.prescribedFlux[key_i][key_j][key_k] = value_k * 0

        if self.realObs:
            self.prescribedConc = self.model(prescribedInit, self.prescribedFlux)
            self.prescribedConc = self.extractObs(self.prescribedConc)
        else:
            self.prescribedConc = self.model(prescribedInit, self.prescribedFlux)

        self.prescribedConcVector = self.o2v(self.prescribedConc)

        return self.prescribedFlux, self.prescribedConcVector

    
    def f2v(self, newInit, fluxes): 
        """
        Returns a vector with the structure of the control vector.
        """
        vector = array(())

        for key_i, value_i in newInit.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in newInit[key_i].items():
                    vector = append(vector, value_j)

        for key_i, value_i in fluxes.items():
            if key_i == 'time' or key_i == 'units':
                pass
            else:
                for key_j, value_j in fluxes[key_i].items():
                    for key_k, value_k in fluxes[key_i][key_j].items():
                        if key_k in self.presFluxKey[key_i]:
                            pass
                        else:
                            vector = append(vector, value_k)
        
        return vector


    def v2f(self, vector): 
        """
        Returns flux and initial condition from a vector.
        """        
        newInit = deepcopy(self.trueInit)

        fluxes = deepcopy(self.priorFluxes)

        i = 0
        for key_i, value_i in newInit.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in newInit[key_i].items():
                    newInit[key_i][key_j] = vector[i]
                    i += 1

        j = 0
        for key_i, value_i in fluxes.items():
            if key_i == 'time' or key_i == 'units':
                pass
            else:
                for key_j, value_j in fluxes[key_i].items():
                    for key_k, value_k in fluxes[key_i][key_j].items():
                        if key_k in self.presFluxKey[key_i]:
                            pass
                        else:
                            fluxes[key_i][key_j][key_k] = vector[i+j*(len(fluxes['time'])):i+(j+1)*(len(fluxes['time']))].reshape(-1,1)
                            j += 1      
        
        return newInit, fluxes


    def genObsVector(self): 
        """
        Returns the true observation vector.
        """
        if self.realObs:
            try:
                self.obsVector = self.o2v(self.trueObs)
            except ValueError:
                if self.trueObs == None:
                    print("An observation set must be enter when 'realObs' is set to True.")
                else:
                    print("Observation set with invalid format.")
        else:
            self.trueObs = self.model(self.trueInit, self.trueFluxes)
            self.obsVector = self.o2v(self.trueObs)

        return self.trueObs, self.obsVector


    def o2v(self, obs):
        """
        Returns a vector with the structure of the observation vector.
        """
        obsVector = array(())

        for key_i, value_i in obs.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in obs[key_i].items():
                    for key_k in obs[key_i][key_j].keys():
                        if key_k == 'time':
                            pass
                        else:
                            obsVector = append(obsVector, obs[key_i][key_j][key_k])

        return obsVector


    def v2o(self, vector):
        """
        Returns observation vector as observations.
        """
        obs = deepcopy(self.trueObs)

        i = 0
        for key_i, value_i in obs.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in obs[key_i].items():
                    obs[key_i][key_j]['value'] = vector[i * len(obs[key_i][key_j]['time']):(i+1) * len(obs[key_i][key_j]['time'])]
                    i += 1

        return obs
    
    
    def joinWithPrescribed(self, fluxes):
        """
        Joins the control vector with the prescribed fluxes.
        """

        completeFlux = deepcopy(self.prescribedFlux)

        for key_i, value_i in fluxes.items():
            if key_i == 'time' or key_i == 'units':
                pass
            else:
                for key_j, value_j in fluxes[key_i].items():
                    for key_k, value_k in fluxes[key_i][key_j].items():
                        if key_k in self.presFluxKey[key_i]:
                            pass
                        else:
                            completeFlux[key_i][key_j][key_k] = value_k

        return completeFlux
    
    
    def genErrorVector(self):
        """
        Generates control and observation error vectors.
        """
        if self.realError:
            self.errorStateVector = self.f2v(self.errorInit, self.errorFluxes)
            self.errorObsVector = self.o2v(self.errorObs)

        else:
            for key_i, value_i in self.trueInit.items():
                if key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.trueInit[key_i].items():
                        self.errorInit[key_i][key_j] = value_j * self.errorInit[key_i][key_j]

            for key_i, value_i in self.priorFluxes.items():
                if key_i == 'time' or key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.priorFluxes[key_i].items():
                        for key_k, value_k in self.priorFluxes[key_i][key_j].items():
                                self.errorFluxes[key_i][key_j][key_k] = value_k * self.errorFluxes[key_i][key_j][key_k]

            self.errorStateVector = self.f2v(self.errorInit, self.errorFluxes)

            for key_i, value_i in self.trueObs.items():
                if key_i == 'units':
                    pass
                else:
                    for key_j, value_j in self.trueObs[key_i].items():
                        for key_k, value_k in self.trueObs[key_i][key_j].items():
                            if key_k == 'time':
                                pass
                            else:
                                self.errorObs[key_i][key_j][key_k] = self.errorObs[key_i][key_j][key_k] * value_k

            self.errorObsVector = self.o2v(self.errorObs)

        return self.errorStateVector, self.errorObsVector


    def extractObs(self, obs):
        """
        Extracts the observations generated by the transport model when the number of observations is inferior than the number of fluxes.
        """

        for key_i, value_i in obs.items():
            if key_i == 'units':
                pass
            else:
                for key_j, value_j in obs[key_i].items():
                    o = array(())
                    obs[key_i][key_j]['time'] = self.trueObs[key_i][key_j]['time']
                    for i in range(len(self.trueFluxes['time'])):
                        if self.trueFluxes['time'][i] in self.trueObs[key_i][key_j]['time']:
                            o = append(o, obs[key_i][key_j]['value'][i])
                            if i == len(self.trueFluxes['time'])-1:
                                obs[key_i][key_j]['value'] = o.reshape(-1,1)
                        else:
                            pass

        return obs
        
        
    def cost(self, posteriorVector):
        """
        Returns the cost function.
        """
        dx = posteriorVector - self.stateVector
        dy = matmul(self.K, posteriorVector) - (self.obsVector - self.prescribedConcVector)
        
        return matmul(matmul(dx.T, inv(self.B)), dx) * 0.5 + matmul(matmul(dy.T, inv(self.R)), dy) * 0.5

    
    def inversion(self):
        """
        Contains the construction of K and G matrices, error covariance matrices B and R,
        and returns the values of the posterior, cost function and matrices.
        """
        self.B = eye(len(self.stateVector))*(self.errorStateVector**2) # Error covariance matrix for state vector
        
        self.R = eye(len(self.obsVector))*(self.errorObsVector**2) # Error covariance matrix for obs
        
        self.K = zeros((len(self.obsVector), len(self.stateVector))) # Kernel matrix
        for i in range(len(self.stateVector)):
            vec = zeros((len(self.stateVector)))
            vec[i] = 1
            c01, flux1 = self.v2f(vec)
            conc1 = self.model(c01, flux1)
            if self.realObs:
                conc1 = self.extractObs(conc1)
            cvec = self.o2v(conc1)
            self.K[:, i] = cvec

        G = matmul(matmul(self.B, self.K.T), inv(matmul(matmul(self.K, self.B), self.K.T) + self.R))
        
        posteriorVector = self.stateVector + matmul(G, ((self.obsVector - self.prescribedConcVector) - matmul(self.K, self.stateVector)))
        posteriorInit, self.posteriorFluxes = self.v2f(posteriorVector)
        completePosterior = self.joinWithPrescribed(self.posteriorFluxes)
        concPosterior = self.model(posteriorInit, completePosterior)
        if self.realObs:
            concPosterior = self.extractObs(concPosterior)
        
        # Cost function values
        costTrue = self.cost(self.trueVector)
        costPrior = self.cost(self.stateVector)
        costPosterior = self.cost(posteriorVector)
        print('Cost function true:', costTrue,'Cost function prior:', costPrior, 'Cost function posterior:', costPosterior)
        
        completePrior = self.joinWithPrescribed(self.priorFluxes)
        concPrior = self.model(self.trueInit, completePrior)
        if self.realObs:
            concPrior = self.extractObs(concPrior)
        
        return completePrior, completePosterior, concPrior, concPosterior