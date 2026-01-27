package aurelia

import spinal.core._
import spinal.lib._

object AureliaTypes {

  // =========================================================
  // SOTA Numeric Types (E8M17 for DSP48 Optimization)
  // =========================================================
  
  // Configuration for Custom Float: 1 Sign, 8 Exp, 17 Mantissa (Implicit 1) = 26 bits
  // This fits perfectly into Xilinx DSP48 (27x18 multiplier)
  def EXP_WIDTH = 8
  def MAN_WIDTH = 17 
  
  // Abstract Numeric Trait for Swapping (Fixed vs Float)
  trait AureliaNumeric[T <: Data] {
    def +(that: T): T
    def -(that: T): T
    def *(that: T): T
    def fromDouble(d: Double): T
    def zero: T
    def one: T
  }

  // Custom Floating Point Implementation
  case class CustomFloat(expWidth: Int, manWidth: Int) extends Bundle with AureliaNumeric[CustomFloat] {
    val sign = Bool()
    val exp  = UInt(expWidth bits)
    val man  = UInt(manWidth bits) // Explicit mantissa (normalized)

    // TODO: Implement actual floating point logic in a separate Math Lib or utilize Spinal's Floating if compatible.
    // For now, we define the structure.
    
    override def +(that: CustomFloat): CustomFloat = {
      // Placeholder: Logic would go here or be inferred by a blackboard/function
      // using a structural implementation
      val res = CustomFloat(expWidth, manWidth)
      res.sign := this.sign ^ that.sign // Dummy Logic for IO connectivity
      res.exp  := this.exp
      res.man  := this.man
      res
    }

    override def -(that: CustomFloat): CustomFloat = {
       val res = CustomFloat(expWidth, manWidth)
       res := this // Dummy
       res
    }

    override def *(that: CustomFloat): CustomFloat = {
       val res = CustomFloat(expWidth, manWidth)
       res := this // Dummy
       res
    }
    
    override def fromDouble(d: Double): CustomFloat = {
       // Constant conversion logic needed here
       val res = CustomFloat(expWidth, manWidth)
       res.sign := False
       res.exp := 127
       res.man := 0
       res
    }
    
    override def zero: CustomFloat = {
      val c = CustomFloat(expWidth, manWidth)
      c.sign := False; c.exp := 0; c.man := 0
      c
    }
    
    override def one: CustomFloat = {
      val c = CustomFloat(expWidth, manWidth)
      c.sign := False; c.exp := (1<<(expWidth-1)).toUInt; c.man := 0
      c
    }
  }

  // Global Type Alias - Switch this to change entire design precision
  // Usage: val a = RealType()
  def RealType() = SFix(16 exp, -16 exp) // Keeping SFix for now until Float logic is fully populated
  // def RealType() = CustomFloat(EXP_WIDTH, MAN_WIDTH) // Uncomment to switch

  def Vector3() = Vec(RealType(), 3)
  def Matrix3() = Vec(Vec(RealType(), 3), 3)
  def Tensor3() = Vec(Matrix3(), 3)
  def Tensor4() = Vec(Tensor3(), 3)
  def StencilWindow() = Vec(Matrix3(), 5)

  def RealConst(value: Double): SFix = { // Update return type if switching
    val ret = RealType()
    ret := value
    ret
  }
  
  def PI_C       = 3.14159265358979323846
  def ONE_VAL    = 1.0
  def ZERO_VAL   = 0.0
  def TWO_VAL    = 2.0
  def EIGHT_VAL  = 8.0
  def TWELVE_VAL = 12.0

  def GRID_DX_C = 1.0e-4
  def TIME_DT_C = 0.01

  def COEFF_ORDER4_1 = 8.0
  def COEFF_ORDER4_2 = -1.0
  
  def INV_12H_C = 1.0 / (TWELVE_VAL * GRID_DX_C)

  def TISSUE_MOBILITY_C     = 0.05
  def MECHANO_SENSITIVITY_C = 1.0e-3

  val NX_C = 20 
  val NY_C = 20
  val NZ_C = 10
  
  def SymmetricMatrix3() = Vec(RealType(), 6)
  
  // Helper for serialization
  def typeWidth = 32 // Approximate for SFix
  
  case class Axi4StreamData(width: Int = 512) extends Bundle {
     val data = Bits(width bits)
     val last = Bool()
     val keep = Bits(width/8 bits)
  }
  
  def Window3D() = Vec(Vec(Vec(Matrix3(), 3), 3), 3)
}
