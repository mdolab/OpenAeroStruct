import numpy as np

import openmdao.api as om


class FailureKS(om.ExplicitComponent):
    """
    Aggregate failure constraints from the structure.

    To simplify the optimization problem, we aggregate the individual
    elemental failure constraints using a Kreisselmeier-Steinhauser (KS)
    function.

    The KS function produces a smoother constraint than using a max() function
    to find the maximum point of failure, which produces a better-posed
    optimization problem.

    The rho inputeter controls how conservatively the KS function aggregates
    the failure constraints. A lower value is more conservative while a greater
    value is more aggressive (closer approximation to the max() function).

    parameters
    ----------
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    Returns
    -------
    failure : float
        KS aggregation quantity obtained by combining the failure criteria
        for each FEM node. Used to simplify the optimization problem by
        reducing the number of constraints.

    """

    def initialize(self):
        self.options.declare("surface", types=dict)
        self.options.declare("rho", types=float, default=100.0)

    def setup(self):
        surface = self.options["surface"]
        rho = self.options["rho"]

        if surface["fem_model_type"] == "tube":
            num_failure_criteria = 2
        elif surface["fem_model_type"] == "wingbox":
            num_failure_criteria = 4
        
        # =================================
        # Adding Tsai Wu here
        # =================================
        elif surface["fem_model_type"] == "tsaiwu_wingbox":
            num_failure_criteria = 16
        # =================================

        self.ny = surface["mesh"].shape[1]

        # =================================
        # Adding an if statement for vonmises
        # =================================
        if surface["fem_model_type"] == "tube" or surface["fem_model_type"] == "wingbox":
            self.add_input("vonmises", val=np.zeros((self.ny - 1, num_failure_criteria)), units="N/m**2")

        # self.add_input("vonmises", val=np.zeros((self.ny - 1, num_failure_criteria)), units="N/m**2")
        # =================================

        # =================================
        # Adding Tsai Wu here
        # =================================
        if surface["fem_model_type"] == "tsaiwu_wingbox":
            self.add_input("tsaiwu_sr", val=np.zeros((self.ny - 1, 16)), units=None)
            self.composite_safetyfactor = surface["composite_safetyfactor"]
            self.srlimit = 1 / self.composite_safetyfactor
        # =================================
        
        self.add_output("failure", val=0.0)
        # self.composite_safetyfactor = surface["composite_safetyfactor"]
        self.sigma = surface["yield"]
        self.rho = rho

        # # =================================
        # # Adding Tsai Wu SF here (srlimit) = 1 / SF
        # # =================================
        # if surface["fem_model_type"] == "tsaiwu_wingbox":
        #     self.srlimit = 1 / surface["composite_safetyfactor"] #NOTE: This option needs to be added in the surface dictionary and connected
        # # =================================

        self.declare_partials("*", "*")

    def compute(self, inputs, outputs):
        sigma = self.sigma
        rho = self.rho
        # =================================
        # Adding Tsai Wu SF here
        # =================================
        # srlimit = 1 / self.composite_safetyfactor
        # =================================

        # =================================
        # Adding an if statement for vonmises calculations and tsaiwu calculations
        # =================================

        # defining the stress_array and stress_limit variables for both cases
        if self.options["surface"]["fem_model_type"] == "tube" or self.options["surface"]["fem_model_type"] == "wingbox":
            stress_array = inputs["vonmises"]
            stress_limit = sigma
        elif self.options["surface"]["fem_model_type"] == "tsaiwu_wingbox":
            stress_array = inputs["tsaiwu_sr"]
            srlimit = 1 / self.composite_safetyfactor
            stress_limit = srlimit

        fmax = np.max(stress_array / stress_limit - 1)

        nlog, nsum, nexp = np.log, np.sum, np.exp
        ks = 1 / rho * nlog(nsum(nexp(rho * (stress_array / stress_limit - 1 - fmax))))
        outputs["failure"] = fmax + ks

        # vonmises = inputs["vonmises"]

        # fmax = np.max(vonmises / sigma - 1)

        # nlog, nsum, nexp = np.log, np.sum, np.exp
        # ks = 1 / rho * nlog(nsum(nexp(rho * (vonmises / sigma - 1 - fmax))))
        # outputs["failure"] = fmax + ks

        # =================================


    def compute_partials(self, inputs, partials):

        # =================================
        # Adding an if statement for vonmises calculations and tsaiwu calculations
        # =================================

        # defining the stress_array and stress_limit variables for both cases
        if self.options["surface"]["fem_model_type"] == "tube" or self.options["surface"]["fem_model_type"] == "wingbox":
            stress_array = inputs["vonmises"]
            stress_limit = self.sigma
        elif self.options["surface"]["fem_model_type"] == "tsaiwu_wingbox":
            stress_array = inputs["tsaiwu_sr"]
            stress_limit = self.srlimit
        
        fmax = np.max(stress_array / stress_limit - 1)
        i, j = np.where((stress_array / stress_limit - 1) == fmax)
        i, j = i[0], j[0]

        ksb = 1.0

        tempb0 = ksb / (self.rho * np.sum(np.exp(self.rho * (stress_array / stress_limit - fmax - 1))))
        tempb = np.exp(self.rho * (stress_array / stress_limit - fmax - 1)) * self.rho * tempb0
        fmaxb = ksb - np.sum(tempb)

        derivs = tempb / stress_limit
        derivs[i, j] += fmaxb / stress_limit

        if self.options["surface"]["fem_model_type"] == "tube" or self.options["surface"]["fem_model_type"] == "wingbox":
            partials["failure", "vonmises"] = derivs.reshape(1, -1)
        elif self.options["surface"]["fem_model_type"] == "tsaiwu_wingbox":
            partials["failure", "tsaiwu_sr"] = derivs.reshape(1, -1)

        # vonmises = inputs["vonmises"]
        # sigma = self.sigma
        # rho = self.rho

        # # Find the location of the max stress constraint
        # fmax = np.max(vonmises / sigma - 1)
        # i, j = np.where((vonmises / sigma - 1) == fmax)
        # i, j = i[0], j[0]

        # # Set incoming seed as 1 so we simply get the jacobian entries
        # ksb = 1.0

        # # Use results from the AD code to compute the jacobian entries
        # tempb0 = ksb / (rho * np.sum(np.exp(rho * (vonmises / sigma - fmax - 1))))
        # tempb = np.exp(rho * (vonmises / sigma - fmax - 1)) * rho * tempb0
        # fmaxb = ksb - np.sum(tempb)

        # # Populate the entries
        # derivs = tempb / sigma
        # derivs[i, j] += fmaxb / sigma

        # # Reshape and save them to the jac dict
        # partials["failure", "vonmises"] = derivs.reshape(1, -1)

        # =================================