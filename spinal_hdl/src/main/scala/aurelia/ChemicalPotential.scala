package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class HighPrecExpUnit extends Component {
  val io = new Bundle {
    val input  = slave Stream(RealType())
    val output = master Stream(RealType())
  }

  val LN2     = RealConst(0.6931471805599453)
  val INV_LN2 = RealConst(1.4426950408889634)
  
  val C2 = RealConst(0.5)
  val C3 = RealConst(0.16666666666666666)
  val C4 = RealConst(0.04166666666666666)
  val C5 = RealConst(0.00833333333333333)

  val s1 = io.input.translateWith {
    val x = io.input.payload
    val k_calc = x * INV_LN2
    
    val k_rounded = k_calc + RealConst(0.5)
    val k_int = k_rounded.asSInt
    
    val k_fix = RealType()
    k_fix := k_int
    
    val r = x - (k_fix * LN2)
    
    (k_int, r)
  }.s2mPipe() 

  val s2 = s1.translateWith {
    val (k, r) = s1.payload
    
    val r2 = r * r
    val r3 = r2 * r
    val r4 = r3 * r
    val r5 = r4 * r
    
    val poly = RealConst(1.0) + r + (r2 * C2) + (r3 * C3) + (r4 * C4) + (r5 * C5)
    
    (k, poly)
  }.s2mPipe()

  val s3 = s2.translateWith {
    val (k, poly) = s2.payload
    
    val result = RealType()
    
    val k_clamped = k.resize(6) 
    
    result := poly |<< k_clamped

    val is_neg = k < 0
    val shift_val = k.abs.asUInt.resize(6)
    
    val shifted = is_neg ? (poly >> shift_val) | (poly << shift_val)
    
    shifted
  }.s2mPipe()

  io.output << s3
}


class ChemicalPotential extends Component {
  val io = new Bundle {
    val ricci_scalar = slave Stream(RealType())
    val mu_out       = master Stream(RealType())
  }

  val MU0_VAL       = RealConst(8.3e-20)
  val SAFE_MU0      = RealConst(1.0e-3)
                                        
  val INV_KT        = RealConst(1.0/4.047e-21)
                                               
  val MAX_SCALE_VAL = RealConst(10.0)
  
  val s1 = io.ricci_scalar.translateWith {
    val r_val = io.ricci_scalar.payload.abs
    val act_val = RealConst(MECHANO_SENSITIVITY_C) * r_val
    
    val scale_val = act_val * RealConst(1.0)
    
    val clamped = (scale_val > MAX_SCALE_VAL) ? MAX_SCALE_VAL | scale_val
    clamped
  }.s2mPipe()
  
  val uExp = new HighPrecExpUnit()
  uExp.io.input << s1
  
  val s3 = uExp.io.output.translateWith {
    val exp_val = uExp.io.output.payload
    val mu = SAFE_MU0 * exp_val
    mu
  }.s2mPipe()
  
  io.mu_out << s3
}
