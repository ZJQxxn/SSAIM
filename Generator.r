'''
Description: Generate multivariate data for experiments of SSAIM.

Programmer: Jiaqi Zhang

Update history: 
    --`2019/04/02: Create the framework of class ``Generator'';
    -- 2019/04/04: Create 7 functions. 
'''
import numpy as np
from FunList import FuncList

class Generator:
    '''
    Description:
        A class for generating simulated datasets for SSAIM.
        
    Attributes:
        -- sampleNum: The number of samples;
        -- varNum: The number of variabels;
        -- s: Sparsity level (number of non-zero elements if theta);
        -- X: Observed sample predictors;
        -- Y: Ovserved responses;
        -- F: Latent variables;
        -- theta: Sparse weight matrix;
        -- funcList: A list of pointers of monotone functions.
        
    Functions:
        (-) __init__: Initialization;
        (-) _zeroFunc: Function correspondings to zero elements of theta;
        (-) _generateTheta: Generate the spars eweight matrix 'theta';
        (-) _generateX: Generate 'X';
        (-) _generateY: Generate 'Y';
        (-) _generateF: Generate 'F';
        (+) generate: Generate the simulated datasets.
    '''

    def __init__(self,sampleNum,varNum,s):
        '''
        Initialization.
        :param sampleNum: The number of samples. 
        :param varNum: The number of variables.
        :param s: Sparsity level.
        '''
        self.sampleNum=sampleNum
        self.varNum=varNum
        self.s=s
        self.X=None
        self.Y=None
        self.theta=np.zeros((self.varNum,1))
        self.F=None
        self.funcList=FuncList().getFunc(self.s)

    def _zeroFunc(self,x):
        '''
        Function corresponds to the indices of zero elements in 'theta'.
        :param x: x.
        :return: Function value, x.
        '''
        return x

    def _generateTheta(self):
        '''
        Generate 'theta'.
        :return: None.
        '''
        for index in range(self.s):
            self.theta[index]=1

    def _generateX(self):
        '''
        Geenrate 'X'.
        :return: None.
        '''
        self.X=np.zeros((self.sampleNum,self.varNum))
        for sampIndex in range(self.sampleNum):
            for varIndex in range(self.varNum):
                if varIndex < self.s:
                    self.X[sampIndex,varIndex]=self.funcList[varIndex](self.F[sampIndex,varIndex])
                else:
                    self.X[sampIndex,varIndex]=self._zeroFunc(self.F[sampIndex,varIndex])

    def _generateY(self):
        '''
        Generate 'Y'.
        :return: None.
        '''
        self.Y=np.zeros(self.sampleNum)
        for sampIndex in range(self.sampleNum):
            epsilon=np.random.multivariate_normal(np.zeros(1),0.25*np.identity(1),1)
            self.Y[sampIndex]=epsilon+self.theta.T @ self.F[sampIndex,:]

    def _generateF(self):
        '''
        Generate 'F'.
        :return: None.
        '''
        # A is generated from a normal distribution
        A=np.random.multivariate_normal(np.zeros(self.varNum),np.identity(self.varNum),self.varNum)
        A_normed = A
        for i in range(A.shape[0]):
            A_normed[i, :] = A[i, :] / np.max(np.abs(A[i, :]))
        self.sigma=A.T @ A
        # Generate F
        self.F=np.random.multivariate_normal(np.zeros(self.varNum),self.sigma,self.sampleNum)

    def generate(self):
        '''
        Generate a simulated dataset.
        :return: None.
        '''
        self._generateTheta()
        self._generateF()
        self._generateX()
        self._generateY()


if __name__ == '__main__':
    # generate data for smoothness validation
    n = 100
    p = 10
    s = 5
    g = Generator(varNum=p, sampleNum=n, s=s)
    g.generate()
    np.savetxt(str.format('\\data_syn_smooth\\theta\\Theta_%dp_%dn.csv', p, eachN), g.theta, delimiter=',')
    np.savetxt(str.format('\\data_syn_smooth\\X\\X_%dp_%dn.csv', p, eachN), g.X, delimiter=',')
    np.savetxt(str.format('\\data_syn_smooth\\Y\\Y_%dp_%dn.csv', p, eachN), g.Y, delimiter=',')
    np.savetxt(str.format('\\data_syn_smooth\\F\\F_%dp_%dn.csv', p, eachN), g.F, delimiter=',')


    # # n = p
    # p = [200,400,600,800,1000,1500,2000]
    # for eachP in p:
    #     print(eachP,"  variables.")
    #     n = eachP
    #     s = int(eachP * 0.01)
    #     g = Generator(varNum=eachP, sampleNum=n, s=s)
    #     g.generate()
    #     np.savetxt(str.format('\\data_0622\\theta\\Theta_%dp_%dn.csv',eachP,n),g.theta,delimiter=',')
    #     np.savetxt(str.format('\\data_0622\\X\\X_%dp_%dn.csv',eachP,n),g.X,delimiter=',')
    #     np.savetxt(str.format('\\data_0622\\Y\\Y_%dp_%dn.csv',eachP,n),g.Y,delimiter=',')
    #     np.savetxt(str.format('\\data_0622\\F\\F_%dp_%dn.csv',eachP,n),g.F,delimiter=',')
    # # Vary n when p=1000
    # p = 1000
    # s=10
    # n=[100,200,400,600,800,2000]
    # for eachN in n:
    #     print(eachN,"  samples.")
    #     g = Generator(varNum=p, sampleNum=eachN, s=s)
    #     g.generate()
    #     np.savetxt(str.format('\\data_0622\\theta\\Theta_%dp_%dn.csv', p, eachN), g.theta, delimiter=',')
    #     np.savetxt(str.format('\\data_0622\\X\\X_%dp_%dn.csv', p, eachN), g.X, delimiter=',')
    #     np.savetxt(str.format('\\data_0622\\Y\\Y_%dp_%dn.csv', p, eachN), g.Y, delimiter=',')
    #     np.savetxt(str.format('\\data_0622\\F\\F_%dp_%dn.csv', p, eachN), g.F, delimiter=',')