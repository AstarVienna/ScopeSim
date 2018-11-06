SimCADO 1.0 Interface Overview
==============================

Interaction Interfaces
----------------------

.. figure:: Design_Interface_Overview.PNG
    :width: 600

    A slimmed down software design interface diagram for SimCADO 1.0

In the above diagram the main interfaces are shown in red. These can
be divided into input and output interfaces, and external and internal
interfaces. The distinction here is that external interfaces refer to
interfaces where parties outside of the SimCADO team are involved, whereas
internal interfaces are those which are within the scope of the SimCADO software
effort:

**External:**

#.  Data files for describing the effects of the optical train,
#.  Parameters from external pieces of software which control and automated
    simulation run,
#.  SimCADO Output which is to be sent to the other data flow teams:
    i.e. Pipeline and Archive,
#.  Output which should be returned to an external source

**Internal:**

#.  The user's local working directory,

    * Input for generating Source objects,
    * Saving output of locally run simulations,
    * Loading simulation setting files

#.  Instrument package folder on the local machine

    * Loading instrument mode description files
    * Loading instrument effect data

#.  The instrument package server, which hosts verified instrument package files


Software Module Interfaces
--------------------------



.. figure:: Design_Interface_Diagram.PNG
    :width: 600

    The software design interface diagram for SimCADO 1.0 including the major
    software components in SimCADO, as well as the third party code (orange
    packages)





Why is this needed ? Is it required for operations (e.g for Observations
preparation, or taking decisions regarding observations, or data reductions) ?
Please detail the use cases involving the usage of the tool by either ESO
operational staff or science users.

"At the time of the PDR, SIMCADO is not a deliverable to ESO. It is therefore
not part of the PDR process. If ESO takes over SIMCADO, it will need a
dedicated review."

"An evaluation phase is needed before accepting or not the proposal to take
ownership.
This phase should evaluate:
    - the needs for the tool as discussed in RIQ YJU5
    - the maintenance costs for ESO:
        The SIMCADO interfaces with other components: pipeline,
        DFS infrastructure, PSF data file, etc..
        The maintenance effort will also depend on the software quality / size /
        design, which will need to be evaluated
If accepted, the tool would need a proper dedicated design/acceptance review
before the MICADO FDR."