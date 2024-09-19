import numpy as np

import openmdao.api as om


class FailureExact(om.ExplicitComponent):
    """
    Output individual failure constraints on each FEM element.

    parameters
    ----------
    vonmises : ny-1 x 2 numpy array
        von Mises stress magnitudes for each FEM element.
    tsaiwu_sr : ny-1 x 4 * numofplies numpy array
        Tsai-Wu strength ratios for each FEM element (ply at each critical element).

    Returns
    -------
    failure : float
        Array of failure conditions. Positive if element has failed. This entity is defined for either
        failure criteria, vonmises or tsaiwu_sr.

    """

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]
        self.useComposite = "useComposite" in self.options["surface"].keys() and self.options["surface"]["useComposite"]
        if self.useComposite:
            ply_angles = surface["ply_angles"]
            numofplies = len(ply_angles)

        if surface["fem_model_type"] == "tube":
            num_failure_criteria = 2

        elif surface["fem_model_type"] == "wingbox":
            if self.useComposite:  # using the Composite wingbox
                num_failure_criteria = 4 * numofplies  # 4 critical elements * number of plies
                self.srlimit = 1 / surface["safety_factor"]
            else:  # using the Isotropic wingbox
                num_failure_criteria = 4

        self.ny = surface["mesh"].shape[1]

        if "safety_factor" in self.options["surface"].keys():
            self.safety_factor = surface["safety_factor"]
        else:
            self.safety_factor = 1

        self.sigma = surface["yield"] / self.safety_factor

        if self.useComposite:  # using the Composite wingbox
            self.add_input("tsaiwu_sr", val=np.zeros((self.ny - 1, num_failure_criteria)), units=None)
        else:  # using the Isotropic structures
            self.add_input("vonmises", val=np.zeros((self.ny - 1, num_failure_criteria)), units="N/m**2")

        self.add_output("failure", val=np.zeros((self.ny - 1, num_failure_criteria)))

        if self.useComposite:  # using the Composite wingbox
            self.declare_partials(
                "failure", "tsaiwu_sr", val=np.eye(((self.ny - 1) * num_failure_criteria)) / self.srlimit
            )
        else:  # using the Isotropic structures
            self.declare_partials(
                "failure", "vonmises", val=np.eye(((self.ny - 1) * num_failure_criteria)) / self.sigma
            )

    def compute(self, inputs, outputs):
        if "vonmises" in inputs:
            outputs["failure"] = inputs["vonmises"] / self.sigma - 1
        elif "tsaiwu_sr" in inputs:
            outputs["failure"] = inputs["tsaiwu_sr"] / self.srlimit - 1
