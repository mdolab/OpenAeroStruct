import openmdao.api as om


class TotalLiftDrag(om.ExplicitComponent):
    """
    Compute the coefficients of lift (CL) and drag (CD) and dimensional lift (L) and drag (D)
    for the entire aircraft, based on the area-weighted sum of individual surfaces' CLs and CDs.

    Parameters
    ----------
    CL : float
        Coefficient of lift (CL) for one lifting surface.
    CD : float
        Coefficient of drag (CD) for one lifting surface.
    S_ref : float
        Surface area for one lifting surface.
    rho : float
        Freestream density.
    v : float
        Freestream velocity.

    Returns
    -------
    CL : float
        Total coefficient of lift (CL) for the entire aircraft.
    CD : float
        Total coefficient of drag (CD) for the entire aircraft.
    L : float
        Total lift force (L) for the entire aircraft.
    D : float
        Total drag force (D) for the entire aircraft.

    """

    def initialize(self):
        self.options.declare("surfaces", types=list)

    def setup(self):
        for surface in self.options["surfaces"]:
            name = surface["name"]
            self.add_input(name + "_CL", val=1.0, tags=["mphys_result"])
            self.add_input(name + "_CD", val=1.0, tags=["mphys_result"])
            self.add_input(name + "_S_ref", val=1.0, units="m**2", tags=["mphys_coupling"])
            self.declare_partials(["CL", "L"], name + "_CL")
            self.declare_partials(["CD", "D"], name + "_CD")
            self.declare_partials(["CL", "L"], name + "_S_ref")
            self.declare_partials(["CD", "D"], name + "_S_ref")

        self.add_input("S_ref_total", val=1.0, units="m**2", tags=["mphys_input"])
        self.add_input("rho", val=1.0, units="kg/m**3", tags=["mphys_input"])
        self.add_input("v", val=1.0, units="m/s", tags=["mphys_input"])
        self.add_output("CL", val=1.0, tags=["mphys_result"])
        self.add_output("CD", val=1.0, tags=["mphys_result"])
        self.add_output("L", val=1.0, units="N", tags=["mphys_result"])
        self.add_output("D", val=1.0, units="N", tags=["mphys_result"])
        self.declare_partials("CL", "S_ref_total")
        self.declare_partials("CD", "S_ref_total")
        self.declare_partials(["L", "D"], ["rho", "v"])

    def compute(self, inputs, outputs):
        # Compute the weighted CL and CD contributions from each surface,
        # weighted by the individual surface areas
        CL = 0.0
        CD = 0.0
        for surface in self.options["surfaces"]:
            name = surface["name"]
            S_ref = inputs[name + "_S_ref"]
            CL += inputs[name + "_CL"] * S_ref
            CD += inputs[name + "_CD"] * S_ref

        # Before normalizing by total area, compute L and D
        outputs["L"] = CL * 0.5 * inputs["rho"] * inputs["v"] ** 2
        outputs["D"] = CD * 0.5 * inputs["rho"] * inputs["v"] ** 2

        # Normalize by total area to get coefficients
        outputs["CL"] = CL / inputs["S_ref_total"]
        outputs["CD"] = CD / inputs["S_ref_total"]

    def compute_partials(self, inputs, partials):
        # Compute the weighted CL and CD contributions from each surface,
        # weighted by the individual surface areas
        CL = 0.0
        CD = 0.0
        for surface in self.options["surfaces"]:
            name = surface["name"]
            S_ref = inputs[name + "_S_ref"]
            CL += inputs[name + "_CL"] * S_ref
            CD += inputs[name + "_CD"] * S_ref

        S_ref_total = inputs["S_ref_total"]

        partials["L", "rho"] = CL * 0.5 * inputs["v"] ** 2
        partials["D", "rho"] = CD * 0.5 * inputs["v"] ** 2

        partials["L", "v"] = CL * inputs["rho"] * inputs["v"]
        partials["D", "v"] = CD * inputs["rho"] * inputs["v"]

        partials["CL", "S_ref_total"] = -CL / S_ref_total**2
        partials["CD", "S_ref_total"] = -CD / S_ref_total**2

        for surface in self.options["surfaces"]:
            name = surface["name"]
            S_ref = inputs[name + "_S_ref"]
            partials["CL", name + "_CL"] = S_ref / S_ref_total
            partials["CD", name + "_CD"] = S_ref / S_ref_total

            partials["CL", name + "_S_ref"] = inputs[name + "_CL"] / S_ref_total
            partials["CD", name + "_S_ref"] = inputs[name + "_CD"] / S_ref_total

            partials["L", name + "_CL"] = 0.5 * inputs["rho"] * inputs["v"] ** 2 * S_ref
            partials["D", name + "_CD"] = 0.5 * inputs["rho"] * inputs["v"] ** 2 * S_ref

            partials["L", name + "_S_ref"] = 0.5 * inputs["rho"] * inputs["v"] ** 2 * inputs[name + "_CL"]
            partials["D", name + "_S_ref"] = 0.5 * inputs["rho"] * inputs["v"] ** 2 * inputs[name + "_CD"]
