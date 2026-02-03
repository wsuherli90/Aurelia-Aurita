package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class MatrixInverse extends Component {
  val io = new Bundle {
    val input  = slave Stream(Matrix3())
    val output = master Stream(Matrix3())
    val det_out = master Stream(RealType()) 
  }

  case class MinorsPayload() extends Bundle {
    val minors = Matrix3()
    val g_in   = Matrix3()
  }
  
  val s1 = io.input.translateWith {
    val g = io.input.payload
    val pay = MinorsPayload()
    
    pay.g_in := g
    
    pay.minors(0)(0) := (g(1)(1) * g(2)(2)) - (g(1)(2) * g(2)(1))
    pay.minors(0)(1) := (g(1)(0) * g(2)(2)) - (g(1)(2) * g(2)(0))
    pay.minors(0)(2) := (g(1)(0) * g(2)(1)) - (g(1)(1) * g(2)(0))
    
    pay.minors(1)(0) := (g(0)(1) * g(2)(2)) - (g(0)(2) * g(2)(1))
    pay.minors(1)(1) := (g(0)(0) * g(2)(2)) - (g(0)(2) * g(2)(0))
    pay.minors(1)(2) := (g(0)(0) * g(2)(1)) - (g(0)(1) * g(2)(0))
    
    pay.minors(2)(0) := (g(0)(1) * g(1)(2)) - (g(0)(2) * g(1)(1))
    pay.minors(2)(1) := (g(0)(0) * g(1)(2)) - (g(0)(2) * g(1)(0))
    pay.minors(2)(2) := (g(0)(0) * g(1)(1)) - (g(0)(1) * g(1)(0))
    
    pay
  }.s2mPipe() 

  case class DetPayload() extends Bundle {
    val det      = RealType()
    val adjugate = Matrix3()
  }

  val s2 = s1.translateWith {
    val m = s1.payload.minors
    val g = s1.payload.g_in
    val p = DetPayload()
    
    p.det := (g(0)(0) * m(0)(0)) - (g(0)(1) * m(0)(1)) + (g(0)(2) * m(0)(2))
    
    p.adjugate(0)(0) :=  m(0)(0) 
    p.adjugate(0)(1) := -m(1)(0) 
    p.adjugate(0)(2) :=  m(2)(0) 
    
    p.adjugate(1)(0) := -m(0)(1) 
    p.adjugate(1)(1) :=  m(1)(1) 
    p.adjugate(1)(2) := -m(2)(1) 
    
    p.adjugate(2)(0) :=  m(0)(2) 
    p.adjugate(2)(1) := -m(1)(2) 
    p.adjugate(2)(2) :=  m(2)(2) 
    
    p
  }.s2mPipe()

  case class InvPayload() extends Bundle {
    val inv_det  = RealType()
    val adjugate = Matrix3()
  }

  
  val rec_unit = new MathUtils.ReciprocalUnit()
  rec_unit.io.input := s2.payload.det
  rec_unit.io.valid_in := s2.valid
  
  val LATENCY = 6 
  val adj_delayed = Delay(s2.payload.adjugate, LATENCY, init = null) 
  
  val s3_synced = Stream(InvPayload())
  s3_synced.valid := rec_unit.io.valid_out
  s3_synced.payload.inv_det  := rec_unit.io.output
  s3_synced.payload.adjugate := adj_delayed
  
  s2.ready := True 
 
  val s4 = s3_synced.translateWith {
    val adj = s3_synced.payload.adjugate
    val inv_d = s3_synced.payload.inv_det
    
    val inv = Matrix3()
    for(i <- 0 until 3; j <- 0 until 3) {
      inv(i)(j) := adj(i)(j) * inv_d
    }
    inv
  }.s2mPipe()

  io.output << s4
  
  val det_delayed = Delay(s2.payload.det, LATENCY, init = null) 
  
  io.det_out.valid   := RegNext(s3_synced.valid) 
  io.det_out.payload := RegNext(det_delayed)      
}
