import numpy as np

import openmdao.api as om


class FailureExact(om.ExplicitComponent):
    """
    Output individual failure constraints on each FEM element.

    Parameters
    ----------
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    Returns
    -------
    failure[ny-1, 2] : numpy array
        Array of failure conditions. Positive if element has failed.

    """

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]

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
        self.sigma = surface["yield"] #NOTE: does this need to be set to null or removed if not being used?
        # =================================
        # Adding Tsai Wu SF here (srlimit) = 1 / SF
        # =================================
        self.srlimit = 1 / surface["composite_safetyfactor"] #NOTE: This option needs to be added in the surface dictionary and connected
        # =================================

        # =================================
        # Adding an if statement for vonmises and tsaiwu_sr
        # =================================
        if surface["fem_model_type"] == "tube" or surface["fem_model_type"] == "wingbox":
            self.add_input("vonmises", val=np.zeros((self.ny - 1, num_failure_criteria)), units="N/m**2")
        elif surface["fem_model_type"] == "tsaiwu_wingbox":
            self.add_input("tsaiwu_sr", val=np.zeros((self.ny - 1, num_failure_criteria)), units=None)

        # self.add_input("vonmises", val=np.zeros((self.ny - 1, num_failure_criteria)), units="N/m**2")
        # =================================

        self.add_output("failure", val=np.zeros((self.ny - 1, num_failure_criteria)))

        # =================================
        # Adding an if statement for vonmises and tsaiwu_sr
        # =================================
        if surface["fem_model_type"] == "tube" or surface["fem_model_type"] == "wingbox":
            self.declare_partials("failure", "vonmises", val=np.eye(((self.ny - 1) * num_failure_criteria)) / self.sigma)
        elif surface["fem_model_type"] == "tsaiwu_wingbox":
            self.declare_partials("failure", "tsaiwu_sr", val=np.eye(((self.ny - 1) * num_failure_criteria)) / self.srlimit)

        # self.declare_partials("failure", "vonmises", val=np.eye(((self.ny - 1) * num_failure_criteria)) / self.sigma)
        # =================================

    # =================================
    # Adding an if statement for vonmises and tsaiwu_sr
    # =================================
    def compute(self, inputs, outputs):
        if "vonmises" in inputs:
            outputs["failure"] = inputs["vonmises"] / self.sigma - 1
        elif "tsaiwu_sr" in inputs:
            outputs["failure"] = inputs["tsaiwu_sr"] / self.srlimit - 1

    # def compute(self, inputs, outputs):
    #     outputs["failure"] = inputs["vonmises"] / self.sigma - 1
    # =================================
