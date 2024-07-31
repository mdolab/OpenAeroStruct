import numpy as np

import openmdao.api as om


class FailureExact(om.ExplicitComponent):
    """
    Output individual failure constraints on each FEM element.

    Parameters
    ----------
    for Isotropic structures:
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    for Composite wingbox:
    tsaiwu_sr[ny-1, 4 * numofplies] : numpy array
        Tsai-Wu strength ratios for each FEM element (ply at each critical element).

    Returns
    -------
    failure[ny-1, 2] : numpy array
        Array of failure conditions. Positive if element has failed. This entity is defined for either
        failure criteria, vonmises or tsaiwu_sr. # TODO: check this

    """

    def initialize(self):
        self.options.declare("surface", types=dict)

    def setup(self):
        surface = self.options["surface"]
        plyangles = surface["plyangles"]
        numofplies = len(plyangles)

        if surface["fem_model_type"] == "tube":
            num_failure_criteria = 2

        elif surface["fem_model_type"] == "wingbox":
            if "useComposite" in surface.keys() and surface["useComposite"]:  # using the Composite wingbox
                num_failure_criteria = 4 * numofplies  # 4 critical elements * number of plies
            else:  # using the Isotropic wingbox
                num_failure_criteria = 4

        self.ny = surface["mesh"].shape[1]
        self.sigma = surface["yield"]

        self.srlimit = 1 / surface["composite_safetyfactor"]

        if "useComposite" in surface.keys() and surface["useComposite"]:  # using the Composite wingbox
            self.add_input("tsaiwu_sr", val=np.zeros((self.ny - 1, num_failure_criteria)), units=None)
        else:  # using the Isotropic structures
            self.add_input("vonmises", val=np.zeros((self.ny - 1, num_failure_criteria)), units="N/m**2")

        self.add_output("failure", val=np.zeros((self.ny - 1, num_failure_criteria)))

        if "useComposite" in surface.keys() and surface["useComposite"]:  # using the Composite wingbox
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
