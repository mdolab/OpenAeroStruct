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
    vonmises : ny-1 x 2 numpy array
        von Mises stress magnitudes for each FEM element.
    tsaiwu_sr : ny-1 x 4 * numofplies numpy array
        Tsai-Wu strength ratios for each FEM element (ply at each critical element).

    Returns
    -------
    failure : float
        KS aggregation quantity obtained by combining the failure criteria
        for each FEM node. Used to simplify the optimization problem by
        reducing the number of constraints. This entity is defined for either
        failure criteria, vonmises or tsaiwu_sr.
    """

    def initialize(self):
        self.options.declare("surface", types=dict)
        self.options.declare("rho", types=float, default=100.0)

    def setup(self):
        surface = self.options["surface"]
        self.rho = self.options["rho"]

        if "safety_factor" in self.options["surface"].keys():
            self.safety_factor = surface["safety_factor"]
        else:
            self.safety_factor = 1

        self.useComposite = "useComposite" in self.options["surface"].keys() and self.options["surface"]["useComposite"]
        if self.useComposite:
            self.num_plies = len(surface["ply_angles"])
            self.input_name = "tsaiwu_sr"
            self.stress_limit = 1 / self.safety_factor
            self.stress_units = None
        else:
            self.input_name = "vonmises"
            self.stress_limit = surface["yield"] / self.safety_factor
            self.stress_units = "N/m**2"

        if surface["fem_model_type"] == "tube":
            num_failure_criteria = 2

        elif surface["fem_model_type"] == "wingbox":
            if self.useComposite:  # using the Composite wingbox
                num_failure_criteria = 4 * self.num_plies  # 4 critical elements * number of plies
            else:  # using the Isotropic wingbox
                num_failure_criteria = 4

        self.ny = surface["mesh"].shape[1]

        self.add_input(self.input_name, val=np.zeros((self.ny - 1, num_failure_criteria)), units=self.stress_units)

        self.add_output("failure", val=0.0)

        self.declare_partials("*", "*")

    def compute(self, inputs, outputs):

        stress_array = inputs[self.input_name]

        fmax = np.max(stress_array / self.stress_limit - 1)

        nlog, nsum, nexp = np.log, np.sum, np.exp
        ks = 1 / self.rho * nlog(nsum(nexp(self.rho * (stress_array / self.stress_limit - 1 - fmax))))
        outputs["failure"] = fmax + ks

    def compute_partials(self, inputs, partials):
        stress_array = inputs[self.input_name]

        fmax = np.max(stress_array / self.stress_limit - 1)
        i, j = np.where((stress_array / self.stress_limit - 1) == fmax)
        i, j = i[0], j[0]

        ksb = 1.0

        tempb0 = ksb / (self.rho * np.sum(np.exp(self.rho * (stress_array / self.stress_limit - fmax - 1))))
        tempb = np.exp(self.rho * (stress_array / self.stress_limit - fmax - 1)) * self.rho * tempb0
        fmaxb = ksb - np.sum(tempb)

        derivs = tempb / self.stress_limit
        derivs[i, j] += fmaxb / self.stress_limit

        if self.useComposite:  # using the Composite wingbox
            partials["failure", "tsaiwu_sr"] = derivs.reshape(1, -1)
        else:  # using the Isotropic structures
            partials["failure", "vonmises"] = derivs.reshape(1, -1)
