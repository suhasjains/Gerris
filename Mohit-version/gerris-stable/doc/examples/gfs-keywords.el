(defvar gfs-abbrevs '(
"Adapt"
"AdaptError"
"AdaptFunction"
"AdaptGradient"
"AdaptStreamlineCurvature"
"AdaptThickness"
"AdaptVorticity"
"Advection"
"AdvectionParams"
"ApproxProjectionParams"
"Axi"
"BcAngle"
"BcDirichlet"
"BcE"
"BcFlather"
"BcNavier"
"BcNeumann"
"BcSubcritical"
"BcTide"
"Boundary"
"BoundaryGradient"
"BoundaryInflowConstant"
"BoundaryMpi"
"BoundaryOutflow"
"BoundaryPeriodic"
"Box"
"BubbleFraction"
"BubbleFractionDt"
"BubbleInteractions"
"CartesianGrid"
"Constant"
"DeferredCompilation"
"Define"
"DerivedVariable"
"Diffusion"
"DischargeElevation"
"Domain"
"DomainProjection"
"DropletToParticle"
"ElectroHydro"
"ElectroHydroAxi"
"Event"
"EventBalance"
"EventFilter"
"EventHarmonic"
"EventList"
"EventScript"
"EventStop"
"EventSum"
"EventSumDirection"
"FeedParticle"
"ForceBuoy"
"ForceCoeff"
"ForceDrag"
"ForceLift"
"Function"
"FunctionConstant"
"FunctionMap"
"FunctionSpatial"
"GEdge"
"GenericMetric"
"GenericSurface"
"GfsAdapt"
"GfsAdaptError"
"GfsAdaptFunction"
"GfsAdaptGradient"
"GfsAdaptStreamlineCurvature"
"GfsAdaptThickness"
"GfsAdaptVorticity"
"GfsAdvection"
"GfsAdvectionParams"
"GfsApproxProjectionParams"
"GfsAxi"
"GfsBcAngle"
"GfsBcDirichlet"
"GfsBcE"
"GfsBcFlather"
"GfsBcNavier"
"GfsBcNeumann"
"GfsBcSubcritical"
"GfsBcTide"
"GfsBoundary"
"GfsBoundaryGradient"
"GfsBoundaryInflowConstant"
"GfsBoundaryMpi"
"GfsBoundaryOutflow"
"GfsBoundaryPeriodic"
"GfsBox"
"GfsBubbleFraction"
"GfsBubbleFractionDt"
"GfsBubbleInteractions"
"GfsCartesianGrid"
"GfsConstant"
"GfsDeferredCompilation"
"GfsDefine"
"GfsDerivedVariable"
"GfsDiffusion"
"GfsDischargeElevation"
"GfsDomain"
"GfsDomainProjection"
"GfsDropletToParticle"
"GfsElectroHydro"
"GfsElectroHydroAxi"
"GfsEvent"
"GfsEventBalance"
"GfsEventFilter"
"GfsEventHarmonic"
"GfsEventList"
"GfsEventScript"
"GfsEventStop"
"GfsEventSum"
"GfsEventSumDirection"
"GfsFeedParticle"
"GfsForceBuoy"
"GfsForceCoeff"
"GfsForceDrag"
"GfsForceLift"
"GfsFunction"
"GfsFunctionConstant"
"GfsFunctionMap"
"GfsFunctionSpatial"
"GfsGEdge"
"GfsGenericMetric"
"GfsGenericSurface"
"GfsGlobal"
"GfsHydrostaticPressure"
"GfsInit"
"GfsInitFaceValues"
"GfsInitFlowConstant"
"GfsInitFraction"
"GfsInitMask"
"GfsInitStokesWave"
"GfsInitVorticity"
"GfsInitWave"
"GfsLayered"
"GfsLayers"
"GfsMap"
"GfsMapFunction"
"GfsMapProjection"
"GfsMapTransform"
"GfsMetric"
"GfsMetricCubed"
"GfsMetricCubed1"
"GfsMetricLaplace"
"GfsMetricLonLat"
"GfsMetricStretch"
"GfsMetricVariable"
"GfsOcean"
"GfsOutput"
"GfsOutputAdaptStats"
"GfsOutputBalance"
"GfsOutputBoundaries"
"GfsOutputCorrelation"
"GfsOutputDiffusionStats"
"GfsOutputDropletSums"
"GfsOutputErrorNorm"
"GfsOutputGRD"
"GfsOutputLocation"
"GfsOutputObject"
"GfsOutputParticle"
"GfsOutputPotentialStats"
"GfsOutputPovrayDF3"
"GfsOutputPPM"
"GfsOutputProgress"
"GfsOutputProjectionStats"
"GfsOutputScalar"
"GfsOutputScalarHistogram"
"GfsOutputScalarMaxima"
"GfsOutputScalarNorm"
"GfsOutputScalarStats"
"GfsOutputScalarSum"
"GfsOutputSimulation"
"GfsOutputSolidForce"
"GfsOutputSolidStats"
"GfsOutputSpectra"
"GfsOutputSquares"
"GfsOutputStreamline"
"GfsOutputTime"
"GfsOutputTiming"
"GfsOutputView"
"GfsParticle"
"GfsParticleForce"
"GfsParticleList"
"GfsParticleToDroplet"
"GfsParticulate"
"GfsParticulateField"
"GfsPdfParticle"
"GfsPhysicalParams"
"GfsPoisson"
"GfsPorous"
"GfsProjectionParams"
"GfsRefine"
"GfsRefineDistance"
"GfsRefineHeight"
"GfsRefineSolid"
"GfsRefineSurface"
"GfsRefineTerrain"
"GfsRemoveDroplets"
"GfsRemovePonds"
"GfsRiver"
"GfsSimulation"
"GfsSimulationMoving"
"GfsSkewSymmetric"
"GfsSolid"
"GfsSolidMoving"
"GfsSource"
"GfsSourceControl"
"GfsSourceControlField"
"GfsSourceCoriolis"
"GfsSourceCulvert"
"GfsSourceDarcy"
"GfsSourceDiffusion"
"GfsSourceDiffusionExplicit"
"GfsSourceElectric"
"GfsSourceFlux"
"GfsSourceFriction"
"GfsSourceGeneric"
"GfsSourceParticulate"
"GfsSourcePipe"
"GfsSourceScalar"
"GfsSourceTension"
"GfsSourceTensionCSS"
"GfsSourceVelocity"
"GfsSourceViscosity"
"GfsSourceViscosityExplicit"
"GfsSpatialSum"
"GfsStoredMetric"
"GfsSurface"
"GfsSurfaceBc"
"GfsSurfaceBc"
"GfsSurfaceGenericBc"
"GfsSurfaceTerrain"
"GfsTerrain"
"GfsTime"
"GfsVariable"
"GfsVariableAge"
"GfsVariableAverage"
"GfsVariableBoolean"
"GfsVariableCurvature"
"GfsVariableDiagonal"
"GfsVariableDistance"
"GfsVariableFiltered"
"GfsVariableFunction"
"GfsVariableLaplacian"
"GfsVariableMetric"
"GfsVariablePoisson"
"GfsVariablePosition"
"GfsVariableResidual"
"GfsVariableStreamFunction"
"GfsVariableTerrain"
"GfsVariableTracer"
"GfsVariableTracerVOF"
"GfsVariableTracerVOFHeight"
"GfsVariableVOFConcentration"
"GfsWave"
"Global"
"HydrostaticPressure"
"Init"
"InitFaceValues"
"InitFlowConstant"
"InitFraction"
"InitMask"
"InitStokesWave"
"InitVorticity"
"InitWave"
"Layered"
"Layers"
"Map"
"MapFunction"
"MapProjection"
"MapTransform"
"Metric"
"MetricCubed"
"MetricCubed1"
"MetricLaplace"
"MetricLonLat"
"MetricStretch"
"MetricVariable"
"Ocean"
"Output"
"OutputAdaptStats"
"OutputBalance"
"OutputBoundaries"
"OutputCorrelation"
"OutputDiffusionStats"
"OutputDropletSums"
"OutputErrorNorm"
"OutputGRD"
"OutputLocation"
"OutputObject"
"OutputParticle"
"OutputPotentialStats"
"OutputPovrayDF3"
"OutputPPM"
"OutputProgress"
"OutputProjectionStats"
"OutputScalar"
"OutputScalarHistogram"
"OutputScalarMaxima"
"OutputScalarNorm"
"OutputScalarStats"
"OutputScalarSum"
"OutputSimulation"
"OutputSolidForce"
"OutputSolidStats"
"OutputSpectra"
"OutputSquares"
"OutputStreamline"
"OutputTime"
"OutputTiming"
"OutputView"
"Particle"
"ParticleForce"
"ParticleList"
"ParticleToDroplet"
"Particulate"
"ParticulateField"
"PdfParticle"
"PhysicalParams"
"Poisson"
"Porous"
"ProjectionParams"
"Refine"
"RefineDistance"
"RefineHeight"
"RefineSolid"
"RefineSurface"
"RefineTerrain"
"RemoveDroplets"
"RemovePonds"
"River"
"Simulation"
"SimulationMoving"
"SkewSymmetric"
"Solid"
"SolidMoving"
"Source"
"SourceControl"
"SourceControlField"
"SourceCoriolis"
"SourceCulvert"
"SourceDarcy"
"SourceDiffusion"
"SourceDiffusionExplicit"
"SourceElectric"
"SourceFlux"
"SourceFriction"
"SourceGeneric"
"SourceParticulate"
"SourcePipe"
"SourceScalar"
"SourceTension"
"SourceTensionCSS"
"SourceVelocity"
"SourceViscosity"
"SourceViscosityExplicit"
"SpatialSum"
"StoredMetric"
"Surface"
"SurfaceBc"
"SurfaceBc"
"SurfaceGenericBc"
"SurfaceTerrain"
"Terrain"
"Time"
"Variable"
"VariableAge"
"VariableAverage"
"VariableBoolean"
"VariableCurvature"
"VariableDiagonal"
"VariableDistance"
"VariableFiltered"
"VariableFunction"
"VariableLaplacian"
"VariableMetric"
"VariablePoisson"
"VariablePosition"
"VariableResidual"
"VariableStreamFunction"
"VariableTerrain"
"VariableTracer"
"VariableTracerVOF"
"VariableTracerVOFHeight"
"VariableVOFConcentration"
"Wave"
)
"Gerris keywords automatically generated by classes.c.")
(defvar gfs-modules '(
"bubbles"
"calavg"
"calculate_cell_average"
"culvert"
"df3"
"electrohydro"
"fft"
"gfsview"
"layered"
"map"
"ode"
"okada"
"particulates"
"porous"
"skewsymmetric"
"stokes"
"terrain"
"tide"
"topics"
)
"Gerris modules automatically generated by modules.c.")
(provide 'gfs-keywords)