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
    for Isotropic structures:
    vonmises[ny-1, 2] : numpy array
        von Mises stress magnitudes for each FEM element.

    for Composite wingbox:
    tsaiwu_sr[ny-1, 4 * numofplies] : numpy array
        Tsai-Wu strength ratios for each FEM element (ply at each critical element).

    Returns
    -------
    failure : float
        KS aggregation quantity obtained by combining the failure criteria
        for each FEM node. Used to simplify the optimization problem by
        reducing the number of constraints. This entity is defined for either
        failure criteria, vonmises or tsaiwu_sr. # TODO: check this

    """

    def initialize(self):
        self.options.declare("surface", types=dict)
        self.options.declare("rho", types=float, default=100.0)

    def setup(self):
        surface = self.options["surface"]
        rho = self.options["rho"]
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

        if "useComposite" in surface.keys() and surface["useComposite"]:
            self.add_input("tsaiwu_sr", val=np.zeros((self.ny - 1, num_failure_criteria)), units=None)
            self.composite_safetyfactor = surface["composite_safetyfactor"]
            self.srlimit = 1 / self.composite_safetyfactor

        else:
            self.add_input("vonmises", val=np.zeros((self.ny - 1, num_failure_criteria)), units="N/m**2")

        self.add_output("failure", val=0.0)
        self.sigma = surface["yield"]
        self.rho = rho

        self.declare_partials("*", "*")

    def compute(self, inputs, outputs):
        sigma = self.sigma
        rho = self.rho

        if (
            "useComposite" in self.options["surface"].keys() and self.options["surface"]["useComposite"]
        ):  # using the Composite wingbox
            stress_array = inputs["tsaiwu_sr"]
            stress_limit = self.srlimit
        else:  # using the Isotropic structures
            stress_array = inputs["vonmises"]
            stress_limit = sigma

        fmax = np.max(stress_array / stress_limit - 1)

        nlog, nsum, nexp = np.log, np.sum, np.exp
        ks = 1 / rho * nlog(nsum(nexp(rho * (stress_array / stress_limit - 1 - fmax))))
        outputs["failure"] = fmax + ks

    def compute_partials(self, inputs, partials):

        if (
            "useComposite" in self.options["surface"].keys() and self.options["surface"]["useComposite"]
        ):  # using the Composite wingbox
            stress_array = inputs["tsaiwu_sr"]
            stress_limit = self.srlimit
        else:  # using the Isotropic structures
            stress_array = inputs["vonmises"]
            stress_limit = self.sigma

        fmax = np.max(stress_array / stress_limit - 1)
        i, j = np.where((stress_array / stress_limit - 1) == fmax)
        i, j = i[0], j[0]

        ksb = 1.0

        tempb0 = ksb / (self.rho * np.sum(np.exp(self.rho * (stress_array / stress_limit - fmax - 1))))
        tempb = np.exp(self.rho * (stress_array / stress_limit - fmax - 1)) * self.rho * tempb0
        fmaxb = ksb - np.sum(tempb)

        derivs = tempb / stress_limit
        derivs[i, j] += fmaxb / stress_limit

        if (
            "useComposite" in self.options["surface"].keys() and self.options["surface"]["useComposite"]
        ):  # using the Composite wingbox
            partials["failure", "tsaiwu_sr"] = derivs.reshape(1, -1)
        else:  # using the Isotropic structures
            partials["failure", "vonmises"] = derivs.reshape(1, -1)
