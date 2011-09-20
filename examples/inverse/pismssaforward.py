import PISM

from siple.gradient.forward import NonlinearForwardProblem
from siple.gradient.nonlinear import InvertNLCG, InvertIGN
from siple.reporting import msg

from linalg_pism import PISMLocalVector
from math import sqrt

class SSAForwardSolver(PISM.ssa.SSASolver):
  def __init__(self,grid,ice,basal,enthalpyconverter,config=None):
    PISM.ssa.SSASolver.__init__(self,grid,ice,basal,enthalpyconverter,config)

  
  def init_vars(self):
    pismVars = PISM.PISMVars()
    for var in [self.surface,self.thickness,self.bed,self.tauc,
                self.enthalpy,self.ice_mask]:
        pismVars.add(var)
        
    # The SSA instance will not keep a reference to pismVars; it only uses it to extract
    # its desired variables.  So it is safe to pass it pismVars and then let pismVars
    # go out of scope at the end of this method.
    self.ssa.init(pismVars)

    if not self.vel_bc is None:
      self.ssa.set_boundary_conditions(self.bc_mask,self.vel_bc)

    self.ssa.setup_vars()

    self.ssa_init = True
    
  def _constructSSA(self):
    return PISM.SSAFEM_Forward(self.grid,self.basal,self.ice,self.enthalpyconverter,self.config)

WIDE_STENCIL = 2

class PISMSSAForwardProblem(NonlinearForwardProblem,PISM.ssa.SSATestCase):
  
  def __init__(self,Mx,My):
    PISM.ssa.SSATestCase.__init__(self,Mx,My)

    self.solver.init_vars()


  def _constructSSA(self):
    return SSAForwardSolver(self.grid,self.ice,self.basal,self.enthalpyconverter,self.config)

  def F(self, x,out=None,guess=None):
    """
    Returns the value of the forward problem at x.

    Nonlinear problems often make use of an initial guess; this can be provided in 'guess'.

    Storage in 'out', if given, is used for the return value.
    """
    # if not guess is None:
    #   self.solver.set_initial_velocity_guess(guess)
    if out is None:
      out = self.rangeVector()
    self.solver.ssa.set_tauc(x.core())
    self.solver.ssa.solveF(out.core())
    return out

  def T(self,d,out=None):
    """
    Returns the value of the linearization, T, of F, at the point x specified previously in linearizeAt, 
    in the direction d.

    Storage in 'out', if given, is used for the return value.
    """
    if out is None:
      out = self.rangeVector()
    self.solver.ssa.solveT(d.core(),out.core())
    return out

  def TStar(self,r,out=None):
    """
    Let T be the linearization of F at the point x (at the point x specified previously in linearizeAt).  
    Its adjoint is T*.  This method returns the value of T* in the direction r.

    Storage in 'out', if given, is used for the return value.
    """
    if out is None:
      out = self.domainVector()
    self.solver.ssa.solveTStar(r.core(),out.core())
    return out

  def linearizeAt(self,x,guess=None):
    """
    Instructs the class that subsequent calls to T and TStar will be conducted for the given value of x.

    Nonlinear problems often make use of an initial guess; this can be provided in 'guess'.
    """
    self.solver.ssa.set_tauc(x.core())

  def evalFandLinearize(self,x,out=None,guess=None):
    """
    Computes the value of F(x) and locks in a linearization.  Sometimes there are efficiencies that
    can be acheived this way.

    Default implementation simply calls F, then linearizeAt.
    """
    if out is None:
      out = self.rangeVector()
    self.linearizeAt(x)
    self.solver.ssa.solveF(out.core())
    return out
  
  def rangeIP(self,a,b):
    """
    Computes the inner product of two vectors in the range.
    """
    return self.solver.ssa.rangeIP(a.core(),b.core())

  def domainIP(self,a,b):
    """
    Computes the inner product of two vectors in the domain.
    """
    return self.solver.ssa.domainIP(a.core(),b.core())

  def rangeVector(self):
    """Constructs a brand new vector from the range vector space"""
    v = PISM.IceModelVec2V()
    v.create(self.grid,"",True,WIDE_STENCIL)
    return PISMLocalVector(v)

  def domainVector(self):
    """Constructs a brand new vector from the domain vector space"""
    v = PISM.IceModelVec2S()
    v.create(self.grid,"",True,WIDE_STENCIL)
    return PISMLocalVector(v)

class InvertSSANLCG(InvertNLCG):
  """
  Inversion of the map

    F: gamma |-> u

  where u is the solution of the PDE

    -Laplacian u + gamma u = f

  F is a map from L^2 to L^2.
  """

  @staticmethod
  def defaultParameters():
    params = InvertNLCG.defaultParameters()
    return params

  def __init__(self,forward_problem,params=None):
    InvertNLCG.__init__(self,params)
    self.forward_problem = forward_problem

  def forwardProblem(self):
    return self.forward_problem

  def stopConditionMet(self,count,x,Fx,y,r):
    """
    Determines if minimization should be halted (based, e.g. on a Morozov discrepancy principle)

    In: count: current iteration count
        x:     point in domain of potential minimizer.
        Fx:    value of nonlinear function at x
        r:     current residual, i.e. y-F(x)    
    """

    J = sqrt(abs(self.forward_problem.rangeIP(r,r)));

    msg('(%d) J=%g goal=%g',count,J,self.Jgoal)

    if( J < self.Jgoal ):
      msg('Stop condition met')
      return True
    return False

  def initialize(self,x,y,deltaLInf):
    """
    Hook called at the start of solve.  This gives the class a chance to massage the input.
    For example, x and y might be dolfin (finite element) Functions; this method should return 
    the associated dolfin GenericVectors.

    The remaining arguments are passed directly from solve, and can be used for determining the
    final stopping criterion.

    Returns dolfin vectors corresponding to the initial value of x and the desired value of y=F(x).    
    """
    xv = PISMLocalVector(x)
    yv = PISMLocalVector(y)

    self.Jgoal = self.params.mu*deltaLInf

    return (xv,yv)

  def finalize(self,x,y):
    """
    Hook called at the end of 'solve'.  Gives the chance to massage the return values.
    """
    tauc = x.core()
    u = y.core()
    return (tauc,u)

  def solve(self,x,y,deltaLInf):
    """
    Solve the ill posed problem F(x)=y where y is know to an L infinity error deltaLInf
    """
    return InvertNLCG.solve(self,x,y,deltaLInf)


class InvertSSAIGN(InvertIGN):
  """
  Inversion of the map

    F: gamma |-> u

  where u is the solution of the PDE

    -Laplacian u + gamma u = f

  F is a map from L^2 to L^2.
  """

  @staticmethod
  def defaultParameters():
    params = InvertIGN.defaultParameters()
    return params

  def __init__(self,forward_problem,params=None):
    InvertIGN.__init__(self,params)
    self.forward_problem = forward_problem

  def forwardProblem(self):
    return self.forward_problem

  def temper_d(self, x,d,y,r):
    dnorm = d.norm('linf');  xnorm = x.norm('linf')
    if dnorm > 2*xnorm:
      msg('wild change predicted by linear step. scaling')
      d.scale(2*xnorm/dnorm)

  def initialize(self,x,y,deltaLInf):
    """
    Hook called at the start of solve.  This gives the class a chance to massage the input.
    For example, x and y might be dolfin (finite element) Functions; this method should return 
    the associated dolfin GenericVectors.

    The remaining arguments are passed directly from solve, and can be used for determining the
    final stopping criterion.

    Returns dolfin vectors corresponding to the initial value of x and the desired value of y=F(x).    
    """
    xv = PISMLocalVector(x)
    yv = PISMLocalVector(y)

    Jgoal = deltaLInf

    return (xv,yv,Jgoal)

  def finalize(self,x,y):
    """
    Hook called at the end of 'solve'.  Gives the chance to massage the return values.
    """

    tauc = x.core()
    u = y.core()
    return (tauc,u)

  def solve(self,x,y,deltaLInf):
    """
    Solve the ill posed problem F(x)=y where y is know to an L infinity error deltaLInf
    """
    return InvertIGN.solve(self,x,y,deltaLInf)