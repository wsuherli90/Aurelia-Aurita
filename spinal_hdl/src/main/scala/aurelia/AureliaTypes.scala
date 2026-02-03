package aurelia

import spinal.core._
import spinal.lib._

object AureliaTypes {

  def EXP_WIDTH = 8
  def MAN_WIDTH = 17 
  
  // Generic numeric trait for mathematical operations
  trait AureliaNumeric[T <: Data] {
    def +(that: T): T
    def -(that: T): T
    def *(that: T): T
    def fromDouble(d: Double): T
    def zero: T
    def one: T
  }

  // Def for default types - eventually to be replaced by Config injection
  def RealType() = SFix(16 exp, -16 exp)
  def Vector3() = Vec(RealType(), 3)
  def Matrix3() = Vec(Vec(RealType(), 3), 3)
  def Tensor3() = Vec(Matrix3(), 3)
  def Tensor4() = Vec(Tensor3(), 3)
  def StencilWindow() = Vec(Matrix3(), 5)

  def SymmetricMatrix3() = Vec(RealType(), 6)

  def RealConst(value: Double): SFix = {
    val ret = RealType()
    ret := value
    ret
  }
  
  // Physical Constants - These can remain as 'Defaults' or Physics Constants
  def PI_C       = 3.14159265358979323846
  def ONE_VAL    = 1.0
  def ZERO_VAL   = 0.0
  def TWO_VAL    = 2.0
  
  // Simulation params that should be in Config, but keeping defaults here for simplified instantiation
  def GRID_DX_DEFAULT = 1.0e-4
  def TIME_DT_DEFAULT = 0.01

  // Configuration for Grid and Hardware Parameters
  // Expert Pattern: Case Class for verification and generation
  case class AureliaConfig(
     nx: Int = 20,
     ny: Int = 20,
     nz: Int = 10,
     maxWidth: Int = 256,
     dataType: HardType[SFix] = HardType(SFix(16 exp, -16 exp)), // Future proofing for generics
     gridDx: Double = GRID_DX_DEFAULT,
     timeDt: Double = TIME_DT_DEFAULT
  )
  
  def Window3D() = Vec(Vec(Vec(Matrix3(), 3), 3), 3)
  def Window5D() = Vec(Vec(Vec(Matrix3(), 5), 5), 5)
  
  // Axi4Stream Definition
  case class Axi4StreamData(width: Int = 512) extends Bundle {
     val data = Bits(width bits)
     val last = Bool()
     val keep = Bits(width/8 bits)
  }

  // Packing/Unpacking utilities
  def unpackSymmetric(sym: Vec[SFix]): Vec[Vec[SFix]] = {
    val m = Matrix3()
    // Diagonal
    m(0)(0) := sym(0)
    m(1)(1) := sym(1)
    m(2)(2) := sym(2)
    // Off-diagonal (Symmetric)
    m(0)(1) := sym(3); m(1)(0) := sym(3)
    m(0)(2) := sym(4); m(2)(0) := sym(4)
    m(1)(2) := sym(5); m(2)(1) := sym(5)
    m
  }

  def packSymmetric(m: Vec[Vec[SFix]]): Vec[SFix] = {
    val s = SymmetricMatrix3()
    s(0) := m(0)(0)
    s(1) := m(1)(1)
    s(2) := m(2)(2)
    s(3) := m(0)(1) 
    s(4) := m(0)(2)
    s(5) := m(1)(2)
    s
  }
}
